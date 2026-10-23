dbglm_moment_terms = function(x, weight) {
  terms = character()
  for (i in seq_along(x)) for (j in seq_len(i)) {
    terms = c(terms, sprintf("SUM((%s)*(%s)*(%s)) AS h_%d_%d", weight, x[i], x[j], i, j))
  }
  terms
}

dbglm_read_information = function(row, p) {
  H = matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(i)) {
    H[i, j] = H[j, i] = row[[sprintf("h_%d_%d", i, j)]]
  }
  H
}

dbglm_basis = function(H, tol) {
  scale = sqrt(pmax(diag(H), 0))
  keep = integer()
  for (j in seq_len(ncol(H))) {
    if (!is.finite(scale[j]) || scale[j] == 0) next
    indices = c(keep, j)
    C = H[indices, indices, drop = FALSE] / outer(scale[indices], scale[indices])
    R = tryCatch(chol(C), error = function(e) NULL)
    if (!is.null(R) && utils::tail(diag(R), 1)^2 > tol) keep = indices
  }
  if (!length(keep)) stop("No identified GLM design columns remain.", call. = FALSE)
  keep
}

dbglm_solve = function(H, b = NULL) {
  scale = sqrt(diag(H))
  if (any(!is.finite(H)) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("GLM working information is invalid; check identification/separation.", call. = FALSE)
  }
  C = H / outer(scale, scale)
  R = tryCatch(chol(C), error = function(e) NULL)
  if (is.null(R)) stop("GLM working information is singular; check separation.", call. = FALSE)
  inverse = chol2inv(R) / outer(scale, scale)
  beta = if (is.null(b)) NULL else as.numeric(backsolve(R, forwardsolve(t(R), b / scale))) / scale
  list(beta = beta, inverse = inverse)
}

dbglm_evaluate = function(conn, table, x, beta, link) {
  cte = dbglm_working_sql(table, x, beta, link)
  terms = dbglm_moment_terms(x, "w")
  terms = c(terms, sprintf("SUM(w*(%s)*z) AS b_%d", x, seq_along(x)),
            sprintf("%s AS s_%d", vapply(sprintf("(%s)*score", x), function(e) dbglm_exact_sum(conn, e), character(1)), seq_along(x)),
            dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
            "SUM(CASE WHEN NOT isfinite(w) OR w <= 0 OR NOT isfinite(eta) THEN 1 ELSE 0 END) AS bad")
  row = dbglm_get_query(conn, paste0(cte, "SELECT ", paste(terms, collapse = ", "), " FROM working"))
  H = dbglm_read_information(row, length(x))
  b = vapply(seq_along(x), function(i) as.numeric(row[[paste0("b_", i)]]), numeric(1))
  score = vapply(seq_along(x), function(i) as.numeric(row[[paste0("s_", i)]]), numeric(1))
  list(H = H, b = b, score = score, score_max = max(abs(score) / sqrt(diag(H))),
       deviance = row$deviance, bad = row$bad, mu_min = row$mu_min, mu_max = row$mu_max)
}

dbglm_irls = function(conn, table, x, beta, family, control, verbose, separation = "error") {
  current = dbglm_evaluate(conn, table, x, beta, family$link)
  initial_information = sqrt(diag(current$H))
  history = current$deviance
  halvings_total = 0L
  for (iter in seq_len(control$glm_max_iter)) {
    candidate = dbglm_solve(current$H, current$b)$beta
    step = 1
    accepted = FALSE
    for (halvings in 0:control$step_max_halvings) {
      next_beta = beta + step * (candidate - beta)
      next_state = dbglm_evaluate(conn, table, x, next_beta, family$link)
      allowance = 64 * .Machine$double.eps * max(1, abs(current$deviance))
      if (next_state$bad == 0 && is.finite(next_state$deviance) &&
          next_state$deviance <= current$deviance + allowance) {
        accepted = TRUE
        break
      }
      step = step / 2
    }
    if (!accepted) stop("GLM step halving failed to find a valid improving step.", call. = FALSE)
    relative = abs(current$deviance - next_state$deviance) / (0.1 + abs(next_state$deviance))
    beta = next_beta
    current = next_state
    history = c(history, current$deviance)
    halvings_total = halvings_total + halvings
    if (isTRUE(verbose)) message("[dbglm] iteration ", iter, ": deviance=",
                                signif(current$deviance, 8), ", score=", signif(current$score_max, 3))
    boundary = dbglm_boundary(current, family$link)
    if (step==1 && relative<=control$glm_tol && boundary &&
        max(abs(current$score)/initial_information)<=control$score_tol) {
      if(separation=="error") stop("GLM reached the numerical probability boundary; a finite solution is not established.",call.=FALSE)
      warning("Returning boundary probabilities only: coefficients are not a finite MLE; inference is unavailable.",call.=FALSE)
      return(c(current,list(beta=beta,iterations=iter,deviance_history=history,
        step_halvings=halvings_total,finite_mle=FALSE,converged=FALSE)))
    }
    if (step == 1 && relative <= control$glm_tol &&
        is.finite(current$score_max) && current$score_max <= control$score_tol) {
      if (dbglm_boundary(current, family$link)) {
        stop("GLM reached the numerical probability boundary; a finite solution is not established.", call. = FALSE)
      }
      return(c(current, list(beta = beta, iterations = iter,
                            deviance_history = history, step_halvings = halvings_total)))
    }
  }
  stop("GLM did not converge within ", control$glm_max_iter, " iterations.", call. = FALSE)
}

dbglm_ssc = function(ssc) {
  defaults = list(K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE,
                  G.adj = TRUE, G.df = "min", t.df = "min")
  if (!is.null(ssc)) {
    if (!is.list(ssc) || is.null(names(ssc)) || anyDuplicated(names(ssc)) ||
        any(!names(ssc) %in% names(defaults))) stop("Invalid `ssc` list.", call. = FALSE)
    defaults[names(ssc)] = ssc
  }
  for (name in c("K.adj", "G.adj", "K.exact")) {
    if (!is.logical(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` adjustment flags must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (name in c("K.fixef", "G.df", "t.df")) {
    if (!is.character(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` policies must be single character strings.", call. = FALSE)
    }
  }
  if (defaults$K.exact || !defaults$K.fixef %in% c("full", "none", "nonnested") ||
      !defaults$G.df %in% c("min", "conventional") || !defaults$t.df %in% c("min", "conventional")) {
    stop("Unsupported `ssc` settings.", call. = FALSE)
  }
  defaults
}

dbglm_covariance = function(conn, table, x, fit, family, kind, ssc, n) {
  inverse = dbglm_solve(fit$H)$inverse
  p = length(x)
  if (ssc$K.adj && n <= p) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  if (kind == "iid") {
    V = inverse
    if (ssc$K.adj) V = V * (n - 1) / (n - p)
    return(list(matrix = V, clusters = NULL))
  }
  cte = dbglm_working_sql(table, x, fit$beta, family$link)
  if (kind == "hc1") {
    row = dbglm_get_query(conn, paste0(cte, "SELECT ",
      paste(dbglm_moment_terms(x, "score*score"), collapse = ", "), " FROM working"))
    G = NULL
  } else {
    scores = sprintf("SUM((%s)*score) AS s_%d", x, seq_along(x))
    terms = dbglm_moment_terms(paste0("s_", seq_along(x)), "1.0")
    row = dbglm_get_query(conn, paste0(cte,
      ", cluster_scores AS (SELECT __dbglm_cluster, ", paste(scores, collapse = ", "),
      " FROM working GROUP BY __dbglm_cluster) SELECT COUNT(*) AS n_clusters, ",
      paste(terms, collapse = ", "), " FROM cluster_scores"))
    G = as.numeric(row$n_clusters)
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  meat = dbglm_read_information(row, p)
  V = inverse %*% meat %*% inverse
  if (ssc$K.adj) V = V * (if (kind == "hc1") n else n - 1) / (n - p)
  if (kind == "cluster" && ssc$G.adj) V = V * G / (G - 1)
  list(matrix = V, clusters = G)
}


build_model_spec = function(formula, conn, table, fe_vars = character(),
                             factor_refs = NULL, reserved_names = character()) {
  design = glm_sql_model_matrix(
    formula = formula,
    conn = conn,
    table = table,
    expand = "all",
    sep = "_x_",
    fe_vars = fe_vars,
    factor_refs = factor_refs
  )
  public_names = design[["coef_names"]]
  if (anyDuplicated(public_names)) {
    duplicates = unique(public_names[duplicated(public_names)])
    stop(
      "Formula terms produce duplicate coefficient names: ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }

  aliases = character(length(public_names))
  reserved = unique(reserved_names)
  for (i in seq_along(public_names)) {
    candidate = sprintf("dbreg_internal_x_%04d__", i)
    while (candidate %in% reserved) candidate = paste0(candidate, "_")
    aliases[[i]] = candidate
    reserved = c(reserved, candidate)
  }
  list(
    rhs = if (inherits(formula, "Formula")) {
      formula(formula, lhs = 0, rhs = 1)
    } else {
      formula
    },
    source_vars = all.vars(if (inherits(formula, "Formula")) {
      formula(formula, lhs = 0, rhs = 1)
    } else {
      formula
    }),
    columns = data.frame(
      alias = aliases,
      name = public_names,
      sql = design[["select_exprs"]],
      term = design[["terms"]],
      stringsAsFactors = FALSE
    ),
    factor_levels = design[["factor_levels"]],
    factor_refs = factor_refs
  )
}

evaluate_model_spec = function(model_spec, newdata) {
  aliases = model_spec[["columns"]][["alias"]]
  if (all(aliases %in% names(newdata))) {
    matrix = Matrix::Matrix(as.matrix(newdata[, aliases, drop = FALSE]), sparse = TRUE)
    colnames(matrix) = aliases
    return(matrix)
  }
  missing_vars = setdiff(model_spec[["source_vars"]], names(newdata))
  if (length(missing_vars)) {
    stop(
      "newdata is missing required regressor columns: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  for (variable in names(model_spec[["factor_levels"]])) {
    training_levels = model_spec[["factor_levels"]][[variable]]
    values = as.character(newdata[[variable]])
    unseen = setdiff(unique(values[!is.na(values)]), as.character(training_levels))
    if (length(unseen)) {
      stop(
        "newdata contains unseen levels for `", variable, "`: ",
        paste(unseen, collapse = ", "),
        call. = FALSE
      )
    }
    newdata[[variable]] = factor(values, levels = as.character(training_levels))
  }

  matrix = Matrix::sparse.model.matrix(model_spec[["rhs"]], data = newdata)
  if ("(Intercept)" %in% colnames(matrix)) {
    matrix = matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  }
  expected = model_spec[["columns"]][["name"]]
  missing_columns = setdiff(expected, colnames(matrix))
  if (length(missing_columns)) {
    stop(
      "Could not reconstruct fitted design columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  matrix = matrix[, expected, drop = FALSE]
  colnames(matrix) = model_spec[["columns"]][["alias"]]
  matrix
}

glm_sql_model_matrix = function(
  formula,
  conn,
  table,
  expand = c("all", "interactions"),
  sep = "_x_",
  fe_vars = character(),
  factor_refs = NULL
) {
  expand = match.arg(expand)
  if (inherits(formula, "Formula")) {
    formula = formula(formula, lhs = 0, rhs = 1)
  }
  tt = terms(formula)
  term_labels = attr(tt, "term.labels")
  if (!length(term_labels)) {
    return(list(select_exprs = character(), col_names = character(),
                coef_names = character(), terms = character(), factor_levels = list()))
  }
  col_info = glm_get_column_info(conn, table, term_labels, factor_refs)
  result = list(
    select_exprs = character(),
    col_names = character(),
    coef_names = character(),
    terms = character()
  )
  for (term in term_labels) {
    expanded = glm_expand_term(term, col_info, expand, sep, all_terms = unique(c(term_labels, fe_vars)), fe_vars = fe_vars)
    result[["select_exprs"]] = c(result[["select_exprs"]], vapply(expanded, `[[`, character(1), "sql"))
    result[["col_names"]] = c(result[["col_names"]], vapply(expanded, `[[`, character(1), "name"))
    result[["coef_names"]] = c(result[["coef_names"]], vapply(expanded, function(x) {
      if (is.null(x[["display"]])) x[["name"]] else x[["display"]]
    }, character(1)))
    result[["terms"]] = c(result[["terms"]], rep(term, length(expanded)))
  }
  result[["factor_levels"]] = col_info[["levels"]]
  result
}

glm_get_column_info = function(conn, table, term_labels, factor_refs = NULL) {
  vars = unique(unlist(lapply(term_labels, function(term) {
    all.vars(parse(text = term, keep.source = FALSE)[[1]])
  })))
  from_sql = if (grepl("^FROM", table, ignore.case = TRUE)) table else paste("FROM", table)
  schema_query = glue("SELECT {paste(vars, collapse = ', ')} {from_sql} LIMIT 0")
  schema = dbGetQuery(conn, schema_query)
  types = list()
  levels = list()
  classes = setNames(vapply(vars, function(v) class(schema[[v]])[1], character(1)), vars)
  factor_vars = vars[
    classes %in% c("character", "factor", "logical") |
      vars %in% names(factor_refs)
  ]
  for (v in vars) types[[v]] = if (v %in% factor_vars) "factor" else "numeric"

  if (detect_backend(conn)[["name"]] == "spark" && length(factor_vars) > 1L) {
    stack_args = unlist(lapply(factor_vars, function(v) {
      c(as.character(DBI::dbQuoteString(conn, v)), sprintf("CAST(%s AS STRING)", v))
    }), use.names = FALSE)
    level_query = paste0(
      "SELECT factor_name, factor_level FROM (SELECT STACK(",
      length(factor_vars), ", ", paste(stack_args, collapse = ", "),
      ") AS (factor_name, factor_level) ", from_sql, ") t ",
      "WHERE factor_level IS NOT NULL GROUP BY factor_name, factor_level ",
      "ORDER BY factor_name, factor_level"
    )
    stacked = dbGetQuery(conn, level_query)
    for (v in factor_vars) {
      raw_levels = stacked$factor_level[stacked$factor_name == v]
      levels[[v]] = switch(
        classes[[v]],
        "integer" = as.integer(raw_levels),
        "numeric" = as.numeric(raw_levels),
        "logical" = as.logical(raw_levels),
        raw_levels
      )
      if (v %in% names(factor_refs)) {
        ref_idx = match(as.character(factor_refs[[v]]), as.character(levels[[v]]))
        if (is.na(ref_idx)) {
          stop(
            "Reference level `", factor_refs[[v]], "` was not found in regressor `", v, "`.",
            call. = FALSE
          )
        }
        levels[[v]] = levels[[v]][c(ref_idx, setdiff(seq_along(levels[[v]]), ref_idx))]
      }
    }
  } else {
    for (v in factor_vars) {
      lvl_query = glue("SELECT DISTINCT {v} {from_sql} WHERE {v} IS NOT NULL ORDER BY {v}")
      levels[[v]] = dbGetQuery(conn, lvl_query)[[1]]
      if (v %in% names(factor_refs)) {
        ref_idx = match(as.character(factor_refs[[v]]), as.character(levels[[v]]))
        if (is.na(ref_idx)) {
          stop(
            "Reference level `", factor_refs[[v]], "` was not found in regressor `", v, "`.",
            call. = FALSE
          )
        }
        levels[[v]] = levels[[v]][c(ref_idx, setdiff(seq_along(levels[[v]]), ref_idx))]
      }
    }
  }
  level_sql = lapply(levels, function(x) {
    vapply(x, function(value) as.character(DBI::dbQuoteLiteral(conn, value)), character(1))
  })
  list(types = types, levels = levels, level_sql = level_sql)
}

parse_power_term = function(term) {
  expr = tryCatch(parse(text = term, keep.source = FALSE)[[1]], error = function(e) NULL)
  if (is.null(expr) || !is.call(expr) || !identical(as.character(expr[[1]]), "I")) {
    return(NULL)
  }
  power = expr[[2]]
  if (!is.call(power) || !identical(as.character(power[[1]]), "^") ||
      !is.symbol(power[[2]]) || !is.numeric(power[[3]]) ||
      !as.numeric(power[[3]]) %in% c(2, 3)) {
    return(NULL)
  }
  variable = as.character(power[[2]])
  exponent = as.integer(power[[3]])
  list(
    variable = variable,
    exponent = exponent,
    sql = sprintf("POWER(%s, %d)", variable, exponent),
    name = sprintf("%s__pow%d", variable, exponent)
  )
}

split_interaction_term = function(term) {
  expr = parse(text = term, keep.source = FALSE)[[1]]
  flatten = function(x) {
    if (is.call(x) && identical(as.character(x[[1]]), ":")) {
      return(c(flatten(x[[2]]), flatten(x[[3]])))
    }
    paste(deparse(x, width.cutoff = 500L), collapse = "")
  }
  flatten(expr)
}

glm_expand_term = function(term, col_info, expand, sep = "_x_", all_terms = character(), fe_vars = character()) {
  vars = split_interaction_term(term)
  is_interaction = length(vars) > 1
  if (!is_interaction && expand == "interactions" && is.null(parse_power_term(term))) {
    return(list(list(sql = term, name = term, display = term)))
  }
  expansions = lapply(vars, function(v) {
    drop_ref = TRUE
    if (is_interaction && identical(col_info[["types"]][[v]], "factor")) {
      other_term = paste(sort(setdiff(vars, v)), collapse = ":")
      canonical_terms = vapply(all_terms, function(t) {
        paste(sort(split_interaction_term(t)), collapse = ":")
      }, character(1))
      drop_ref = v %in% fe_vars || other_term %in% canonical_terms
    }
    glm_expand_variable(v, col_info, expand, is_interaction, drop_ref)
  })
  glm_cross_product(expansions, sep)
}

glm_expand_variable = function(var, col_info, expand, in_interaction, drop_ref = TRUE) {
  power = parse_power_term(var)
  source_var = if (is.null(power)) var else power[["variable"]]
  is_factor = is.null(power) && identical(col_info[["types"]][[source_var]], "factor")
  if (is_factor && (expand == "all" || in_interaction)) {
    lvls = col_info[["levels"]][[source_var]]
    lvl_sql = col_info[["level_sql"]][[source_var]]
    if (length(lvls) < 2) {
      return(list(list(
        sql = "1", name = paste0(var, "_constant"),
        display = paste0(var, "_constant")
      )))
    }
    if (drop_ref) {
      lvls = lvls[-1]
      lvl_sql = lvl_sql[-1]
    }
    lapply(seq_along(lvls), function(i) {
      lvl = lvls[[i]]
      list(
        sql = glue("CASE WHEN {source_var} = {lvl_sql[[i]]} THEN 1.0 ELSE 0.0 END"),
        name = paste0(source_var, lvl),
        display = paste0(source_var, lvl)
      )
    })
  } else if (!is.null(power)) {
    list(list(
      sql = power[["sql"]],
      name = power[["name"]],
      display = sprintf("I(%s^%d)", power[["variable"]], power[["exponent"]])
    ))
  } else {
    list(list(sql = var, name = var, display = var))
  }
}

glm_cross_product = function(expansions, sep = "_x_") {
  if (length(expansions) == 1) return(expansions[[1]])
  rest = glm_cross_product(expansions[-1], sep)
  result = list()
  for (e1 in expansions[[1]]) {
    for (e2 in rest) {
      e1_sql = e1[["sql"]]
      e2_sql = e2[["sql"]]
      result = c(result, list(list(
        sql = glue("({e1_sql}) * ({e2_sql})"),
        name = paste0(e1[["name"]], sep, e2[["name"]]),
        display = paste0(
          if (is.null(e1[["display"]])) e1[["name"]] else e1[["display"]],
          ":",
          if (is.null(e2[["display"]])) e2[["name"]] else e2[["display"]]
        )
      )))
    }
  }
  result
}
