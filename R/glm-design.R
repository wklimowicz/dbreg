# A single GLM design specification shared by SQL fitting and R prediction.
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
    alias_data = model_matrix_data(~ ., newdata[, aliases, drop = FALSE])
    matrix = Matrix::Matrix(as.matrix(alias_data), sparse = TRUE)
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

  design_data = model_matrix_data(model_spec[["rhs"]], newdata)
  matrix = Matrix::sparse.model.matrix(model_spec[["rhs"]], data = design_data)
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
      raw_levels = stacked[["factor_level"]][stacked[["factor_name"]] == v]
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
