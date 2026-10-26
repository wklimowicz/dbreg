# Small-sample policies count varying-slope coefficients separately from FE intercepts.
fe_intercept_df = function(sizes) {
  if (length(sizes)) sum(sizes) - length(sizes) + 1L else 0L
}

# Variance-covariance and meat matrix computation
# These functions are shared across all dbreg strategies.

#' Compute variance-covariance matrix
#' @keywords internal
compute_vcov = function(
  vcov_type = "iid",
  strategy = "compress",
  XtX_inv,
  rss,
  df_res,
  nobs_orig, # N
  n_params = NULL, # K (for CR1 correction)
  X = NULL,
  rss_g = NULL,
  meat = NULL,
  sandwich_bread = NULL,
  ssc = NULL
) {
  ssc = dbglm_ssc(if (is.null(ssc)) "full" else ssc)
  if (is.null(n_params)) n_params = nobs_orig - df_res
  if (!is.finite(df_res) || df_res <= 0) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  if (!is.finite(rss) || rss < 0) stop("Residual sum of squares is invalid; numerical precision was lost.", call. = FALSE)
  weight_scale = attr(meat, "weight_scale")
  if (is.null(weight_scale)) weight_scale = 1
  if (is.null(sandwich_bread)) sandwich_bread = XtX_inv
  sandwich_bread = sandwich_bread * weight_scale
  if (vcov_type == "hc1") {
    if (strategy == "compress" && is.null(meat)) {
      # Compress strategy: HC1 with grouped residuals
      meat = crossprod(X, Diagonal(x = as.numeric(rss_g)) %*% X)
    }
    # meat should be provided for demean/mundlak/moments strategies
    if (is.null(meat)) {
      stop("HC1 requires meat matrix for non-compress strategies")
    }
    scale_hc1 = if (ssc[["K.adj"]]) nobs_orig / (nobs_orig - n_params) else 1
    vcov_mat = scale_hc1 * (sandwich_bread %*% meat %*% t(sandwich_bread))
    attr(vcov_mat, "type") = "hc1"
  } else if (vcov_type == "cluster") {
    # Cluster-robust (CR1) standard errors
    if (is.null(meat)) {
      stop("Clustered SEs require a meat matrix from compute_score_meat()")
    }
    n_clusters = attr(meat, "n_clusters") # G
    if (is.null(n_clusters)) {
      stop("Meat matrix missing n_clusters attribute")
    }
    if (any(n_clusters < 2)) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    if (is.null(n_params)) n_params = ncol(XtX_inv)
    # CR1 small-sample correction: (G/(G-1)) * ((N-1)/(N-K)), as in fixest
    scale_cr1 = (if (ssc[["G.adj"]]) n_clusters / (n_clusters - 1) else 1) *
      (if (ssc[["K.adj"]]) (nobs_orig - 1) / (nobs_orig - n_params) else 1)
    vcov_mat = scale_cr1 * (sandwich_bread %*% meat %*% t(sandwich_bread))
    attr(vcov_mat, "type") = "cluster"
    attr(vcov_mat, "n_clusters") = n_clusters
  } else {
    # IID case (same for all strategies)
    sigma2 = rss / (if (ssc[["K.adj"]]) nobs_orig - n_params else nobs_orig - 1)
    vcov_mat = sigma2 * XtX_inv
    attr(vcov_mat, "type") = "iid"
  }
  attr(vcov_mat, "n_params") = attr(vcov_mat, "df.K") = n_params
  attr(vcov_mat, "ssc") = ssc
  attr(vcov_mat, "df.t") = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else if (vcov_type == "cluster" && ssc[["t.df"]] == "min") min(n_clusters) - 1 else max(nobs_orig - n_params, 1)
  dimnames(vcov_mat) = dimnames(XtX_inv)
  vcov_mat
}

#' Linear predictor as a SQL expression
#'
#' Pushes fitted coefficients back into the database as `b0 + b1*x1 + ...`.
#' This is the shared prediction scaffolding behind the score pass (residuals
#' for robust/clustered meat matrices), and can equally power out-of-sample
#' predictions or residualization steps in other estimators.
#'
#' @param x_sql Character vector of SQL column expressions.
#' @param beta Numeric coefficients aligned with `x_sql`.
#' @param intercept Optional numeric intercept.
#' @keywords internal
sql_linear_predictor = function(x_sql, beta, intercept = NULL) {
  terms = c(
    if (!is.null(intercept)) sprintf("%.17g", intercept),
    if (length(x_sql)) sprintf("%.17g * %s", beta, x_sql)
  )
  if (!length(terms)) {
    return("0")
  }
  paste0("(", paste(terms, collapse = " + "), ")")
}

#' Compute HC1 or cluster-robust meat matrices via a SQL score pass
#'
#' Aggregates the score vectors `w * e * s` for one or more residuals `e`
#' (SQL expressions, typically `y - x'b` built with [sql_linear_predictor()])
#' that share the same score columns `s`:
#' - `vcov_type = "hc1"`: `sum_i w_i^2 e_i^2 s_i s_i'`
#' - `vcov_type = "cluster"`: sum over clusters g of `s_g s_g'`, where `s_g`
#'   is the sum of `w_i e_i s_i` over observations i in cluster g
#'
#' All residuals are scored in a single query. For OLS the score columns are
#' the regressors; for IV they are the instruments (shared by the 2SLS fit and
#' its first stages), while each residual is formed from its own regressors.
#'
#' @param inputs Inputs environment (for `conn`, `vcov_type_req` and
#'   `cluster_var`). Returns `NULL` when `vcov_type_req` is `"iid"`, since no
#'   meat matrix is needed.
#' @param frame Estimation frame (for `cte_sql`, `cte_name` and
#'   `weights_expr`); see [build_base_frame()].
#' @param resid_sql Named character vector of SQL residual expressions.
#' @param score_sql,score_names SQL expressions and names for the score columns.
#' @param score_intercept Logical. Include a constant score column (first).
#' @return Named list of meat matrices (one per residual) with dimnames
#'   `c("(Intercept)"?, score_names)`. Clustered meats carry `n_clusters`.
#' @keywords internal
compute_score_meat = function(
  inputs,
  frame,
  resid_sql,
  score_sql,
  score_names = score_sql,
  score_intercept = FALSE
) {
  env2env(inputs, environment(), c("conn", "cluster_var"))
  cte_sql = frame[["cte_sql"]]
  cte_name = frame[["cte_name"]]
  weights_expr = frame[["weights_expr"]]
  weight_scale = frame[["weight_scale"]]
  if (is.null(weight_scale)) weight_scale = 1
  if (!is.null(weights_expr)) weights_expr = paste0("(", weights_expr, ")/(", dbglm_number(weight_scale), ")")
  vcov_type = inputs[["vcov_type_req"]]
  if (!vcov_type %in% c("hc1", "cluster")) {
    return(NULL)
  }
  as_double = function(x) sql_cast_double(x, conn)
  resid = as_double(resid_sql)
  s_sql = c(if (score_intercept) "1", score_sql)
  s_names = c(if (score_intercept) "(Intercept)", score_names)
  s_term = ifelse(s_sql == "1", NA_character_, as_double(s_sql))
  idx = seq_along(s_sql)
  sql_prod = function(...) {
    x = c(...)
    paste(x[!is.na(x)], collapse = " * ")
  }

  pairs = which(upper.tri(diag(length(idx)), diag = TRUE), arr.ind = TRUE)
  if (vcov_type == "hc1") {
    w_sq = if (is.null(weights_expr)) NULL else as_double(sprintf("(%s) * (%s)", weights_expr, weights_expr))
    terms = unlist(lapply(seq_along(resid), function(r) {
      sprintf(
        "SUM(%s) AS meat_%d_%d_%d",
        vapply(seq_len(nrow(pairs)), function(k) {
          sql_prod(resid[r], resid[r], s_term[pairs[k, 1]], s_term[pairs[k, 2]], w_sq)
        }, character(1)),
        r, pairs[, 1], pairs[, 2]
      )
    }))
    sql = paste0(cte_sql, ",\nmeat AS (SELECT ", paste(terms, collapse = ", "), " FROM ", cte_name, ")\nSELECT * FROM meat")
  } else {
    w = if (is.null(weights_expr)) NULL else as_double(weights_expr)
    terms = unlist(lapply(seq_along(resid), function(r) {
      sprintf(
        "SUM(%s) AS score_%d_%d",
        vapply(idx, function(i) sql_prod(resid[r], s_term[i], w), character(1)),
        r, idx
      )
    }))
    cross_terms = unlist(lapply(seq_along(resid), function(r) sprintf(
      "SUM(score_%d_%d * score_%d_%d) AS meat_%d_%d_%d",
      r, pairs[, 1], r, pairs[, 2], r, pairs[, 1], pairs[, 2]
    )))
    sql = paste0(
      cte_sql,
      ",\ncluster_scores AS (\n  SELECT ", cluster_var, ", ", paste(terms, collapse = ", "),
      "\n  FROM ", cte_name, "\n  GROUP BY ", cluster_var, "\n)\nSELECT COUNT(*) AS n_clusters, ",
      paste(cross_terms, collapse = ","), " FROM cluster_scores"
    )
  }
  res = dbGetQuery(conn, sql)

  meats = lapply(seq_along(resid), function(r) {
    meat = matrix(0, length(idx), length(idx), dimnames = list(s_names, s_names))
    for (k in seq_len(nrow(pairs))) {
      i = pairs[k, 1]
      j = pairs[k, 2]
      meat[i, j] = meat[j, i] = res[[sprintf("meat_%d_%d_%d", r, i, j)]]
    }
    if (vcov_type == "cluster") attr(meat, "n_clusters") = as.numeric(res[["n_clusters"]])
    attr(meat, "weight_scale") = weight_scale
    meat
  })
  names(meats) = names(resid_sql)
  meats
}

inference_ssc = function(ssc) dbglm_ssc(ssc)

inference_cluster_terms = function(value) {
  if (is.character(value)) {
    if (!length(value) || anyNA(value) || any(!nzchar(value))) stop("Cluster columns must be nonempty names.", call. = FALSE)
    return(lapply(value, function(name) name))
  }
  if (!inherits(value, "formula")) stop("Clusters require column names or a formula.", call. = FALSE)
  if (length(value) == 3L && (!is.symbol(value[[2L]]) || !as.character(value[[2L]]) %in% c("cluster", "twoway"))) {
    stop("Unsupported covariance formula type.", call. = FALSE)
  }
  node = value[[length(value)]]
  split = function(node) {
    if (is.call(node) && identical(node[[1L]], as.name("+"))) return(c(split(node[[2L]]), split(node[[3L]])))
    list(node)
  }
  components = function(node) {
    if (is.symbol(node)) return(as.character(node))
    if (is.call(node) && identical(node[[1L]], as.name("^"))) return(unique(c(components(node[[2L]]), components(node[[3L]]))))
    stop("Cluster terms must be columns or typed column combinations with ^.", call. = FALSE)
  }
  lapply(split(node), components)
}

inference_vcov_request = function(object, vcov, cluster = NULL) {
  if (!is.null(cluster) && !is.null(vcov)) stop("Specify either vcov or cluster.", call. = FALSE)
  if (!is.null(cluster)) return(list(kind = "cluster", clusters = inference_cluster_terms(cluster)))
  if (inherits(vcov, "formula")) return(list(kind = "cluster", clusters = inference_cluster_terms(vcov)))
  if (!is.character(vcov) || length(vcov) != 1L || is.na(vcov)) stop("Unsupported covariance request.", call. = FALSE)
  kind = tolower(vcov)
  if (kind %in% c("iid", "standard")) return(list(kind = "iid", clusters = list()))
  if (kind %in% c("hc1", "hetero", "white")) return(list(kind = "hc1", clusters = list()))
  if (!kind %in% c("cluster", "twoway")) stop("Unsupported covariance type: ", vcov, ".", call. = FALSE)
  groups = object[["fitted_state"]][["inference"]][["fe_components"]]
  count = if (kind == "cluster") 1L else 2L
  if (length(groups) < count) stop("This covariance type requires explicit cluster columns.", call. = FALSE)
  list(kind = "cluster", clusters = groups[seq_len(count)])
}

inference_fe_nesting = function(conn, state, clusters) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  nested = rep(FALSE, length(state[["fe_components"]]))
  if (!length(nested)) return(nested)
  labels = vapply(state[["fe_components"]], paste, collapse = "^", character(1))
  cluster_labels = vapply(clusters, paste, collapse = "^", character(1))
  nested = vapply(seq_along(nested), function(i) labels[[i]] %in% cluster_labels ||
    any(state[["fe_components"]][[i]] %in% cluster_labels), logical(1))
  # Match fixest's K.exact=FALSE shortcut when every cluster is itself an FE.
  # Otherwise check the remaining functional dependencies on the fitted rows.
  if (all(nested) || all(cluster_labels %in% labels)) return(nested)
  queries = unlist(lapply(which(!nested), function(i) {
    vapply(clusters, function(cluster) paste0(
      "SELECT ", i, " AS fe,MAX(n) AS maximum FROM (SELECT ",
      sql_distinct_tuple_count(conn, cluster), " AS n FROM ", q(state[["table"]]), " WHERE ", q(state[["estimated"]]),
      " GROUP BY ", paste(q(state[["fe_components"]][[i]]), collapse = ","), ") g"
    ), character(1))
  }), use.names = FALSE)
  rows = backend_get_query(conn, paste(queries, collapse = " UNION ALL "))
  nested[as.integer(rows[["fe"]][as.numeric(rows[["maximum"]]) == 1])] = TRUE
  nested
}

inference_parameter_count = function(conn, state, request, ssc) {
  K = state[["K_model"]]
  sizes = state[["fe_sizes"]]
  fe_df = state[["fe_df"]]
  if (is.null(fe_df)) fe_df = if (length(sizes)) sum(sizes) - length(sizes) + 1 else 0
  if (ssc[["K.fixef"]] == "none") return(K - fe_df)
  if (request[["kind"]] == "cluster" && ssc[["K.fixef"]] == "nonnested" && length(sizes)) {
    nested = inference_fe_nesting(conn, state, request[["clusters"]])
    K = K - sum(sizes[nested] - 1)
    K = max(K, length(state[["coefficient_names"]]) + 1L)
  }
  K
}

inference_repair_vcov = function(V, fix = TRUE) {
  if (!length(V)) return(V)
  if (any(!is.finite(V))) stop("Post-fit covariance is nonfinite.", call. = FALSE)
  e = eigen(V, symmetric = TRUE)
  diagnostics = list(
    min_eigenvalue = min(e[["values"]]), requested_fix = fix,
    repaired = FALSE, max_adjustment = 0
  )
  if (all(e[["values"]] > 0)) {
    attr(V, "repair") = diagnostics
    return(V)
  }
  repaired = tcrossprod(e[["vectors"]] %*% diag(pmax(e[["values"]], 1e-16), nrow(V)), e[["vectors"]])
  dimnames(repaired) = dimnames(V)
  if (max(abs(V - repaired)) > 1e-8) warning(if (fix)
    "The VCOV matrix is not positive definite and was 'fixed'." else "The VCOV matrix is not positive definite.", call. = FALSE)
  diagnostics[["max_adjustment"]] = max(abs(V - repaired))
  diagnostics[["repaired"]] = fix
  result = if (fix) repaired else V
  attr(result, "repair") = diagnostics
  result
}

retained_covariance = function(object, vcov = "iid", cluster = NULL, ssc = NULL, vcov_fix = TRUE) {
  state = object[["fitted_state"]][["inference"]]
  if (is.null(state)) stop("Post-fit covariance requires fitting with retain='all' or retain='inference'.", call. = FALSE)
  if (!identical(state[["version"]], 1L)) stop("Unsupported retained-inference version.", call. = FALSE)
  if (!is.logical(vcov_fix) || length(vcov_fix) != 1L || is.na(vcov_fix)) stop("vcov_fix must be TRUE or FALSE.", call. = FALSE)
  conn = fitted_resource_connection(object[["resource_state"]])
  request = inference_vcov_request(object, if (is.null(cluster)) vcov else NULL, cluster)
  clusters = request[["clusters"]]
  if (length(clusters) > 2L || (inherits(object, "dbglm") && !isTRUE(object[["linear_model"]]) && length(clusters) > 1L)) {
    stop("Supported clustering is one-/two-way OLS or one-way GLM.", call. = FALSE)
  }
  if (is.null(ssc)) ssc = state[["default_ssc"]]
  ssc = inference_ssc(ssc)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  from = paste("FROM", q(state[["table"]]), "WHERE", q(state[["estimated"]]))
  if (length(clusters)) {
    columns = unique(unlist(clusters))
    classes = state[["column_classes"]]
    if (any(!columns %in% names(classes))) stop("Cluster columns are absent from the retained source.", call. = FALSE)
    missing = paste0(q(columns), " IS NULL")
    numeric = columns[classes[columns] %in% c("numeric", "integer")]
    if (length(numeric)) missing = c(missing, paste0("isnan(CAST(", q(numeric), " AS DOUBLE))"))
    bad = backend_get_query(conn, paste("SELECT COUNT(*) AS n", from, "AND (", paste(missing, collapse = " OR "), ")"))
    if (as.numeric(bad[["n"]]) > 0) stop("Post-fit clusters have missing values on the fitted sample; the sample cannot change.", call. = FALSE)
  }
  n = state[["nobs"]]
  K = inference_parameter_count(conn, state, request, ssc)
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  bread = state[["bread"]]
  score_projection = state[["score_projection"]]
  if (is.null(score_projection)) score_projection = bread
  G = numeric()
  if (request[["kind"]] == "iid") {
    if (state[["kind"]] == "ols") V = bread * (state[["rss"]] / state[["weight_scale"]]) / (if (ssc[["K.adj"]]) n - K else n - 1) else if (state[["kind"]] == "gaussian")
      V = bread * state[["dispersion_scaled"]] * (if (ssc[["K.adj"]]) (n - 1) / (n - K) else 1) else
      V = bread / state[["weight_scale"]] * (if (ssc[["K.adj"]]) (n - 1) / (n - K) else 1)
  } else if (request[["kind"]] == "hc1") {
    V = score_moments(conn, from, q(state[["scores"]]), projection = score_projection)[["matrix"]]
    if (ssc[["K.adj"]]) V = V * n / (n - K)
  } else {
    components = c(clusters, if (length(clusters) == 2L) list(unique(unlist(clusters))))
    moments = lapply(components, function(columns) score_moments(conn, from, q(state[["scores"]]), q(columns), projection = score_projection))
    G = vapply(moments, `[[`, numeric(1), "groups")
    if (any(G < 2)) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    meat = bread * 0
    for (i in seq_along(moments)) {
      adjustment = if (!ssc[["G.adj"]]) 1 else if (ssc[["G.df"]] == "min") min(G[seq_along(clusters)]) / (min(G[seq_along(clusters)]) - 1) else G[[i]] / (G[[i]] - 1)
      meat = meat + (if (i > length(clusters)) -1 else 1) * adjustment * moments[[i]][["matrix"]]
    }
    V = meat
    if (ssc[["K.adj"]]) V = V * (n - 1) / (n - K)
    G = G[seq_along(clusters)]
  }
  V = inference_repair_vcov((V + t(V)) / 2, vcov_fix)
  dimnames(V) = rep(list(state[["coefficient_names"]]), 2L)
  attr(V, "type") = attr(V, "vcov_type") = request[["kind"]]
  attr(V, "cluster_var") = if (length(clusters)) vapply(clusters, paste, collapse = "^", character(1))
  attr(V, "n_params") = attr(V, "df.K") = K
  attr(V, "ssc") = ssc
  attr(V, "rss") = state[["rss"]]
  if (!is.null(attr(object[["vcov"]], "tss"))) attr(V, "tss") = attr(object[["vcov"]], "tss")
  attr(V, "df.t") = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else if (length(G) && ssc[["t.df"]] == "min") max(min(G) - 1, 1) else max(n - K, 1)
  attr(V, "df_t") = attr(V, "df.t")
  if (length(G)) attr(V, "n_clusters") = G
  V
}

# Aggregate cluster score cross-products in SQL. Even observation-sized cluster
# dictionaries return one matrix row; intersections are typed GROUP BY tuples.
sql_score_moments = function(from_statement, score_sql, cluster_sql = character(), projection = NULL) {
  columns = if (length(score_sql)) paste0("dbreg_score_", seq_along(score_sql)) else character()
  values = if (length(columns)) paste(
    if (length(cluster_sql)) paste0("SUM(", score_sql, ")") else score_sql,
    "AS", columns
  ) else "1 AS dbreg_score_empty__"
  source = paste0(
    "(SELECT ", paste(values, collapse = ","), " ", from_statement,
    if (length(cluster_sql)) paste(" GROUP BY", paste(cluster_sql, collapse = ",")) else "", ") dbreg_scores"
  )
  if (!is.null(projection) && length(columns)) {
    stopifnot(nrow(projection) == length(columns), ncol(projection) == length(columns))
    expressions = vapply(seq_len(ncol(projection)), function(i)
      sql_linear_predictor(columns, projection[, i]), character(1))
    projected = paste0("dbreg_influence_", seq_along(columns))
    source = paste0("(SELECT ", paste(expressions, "AS", projected, collapse = ","), " FROM ", source, ") dbreg_influences")
    columns = projected
  }
  paste("SELECT", paste(c("COUNT(*) AS n_groups", dbglm_moment_terms(columns, "1.0")), collapse = ","), "FROM", source)
}

score_moments = function(conn, from_statement, score_sql, cluster_sql = character(), projection = NULL) {
  row = backend_get_query(conn, sql_score_moments(from_statement, score_sql, cluster_sql, projection))
  if (nrow(row) != 1L) stop("Score moment aggregation must return one row.", call. = FALSE)
  list(matrix = dbglm_read_information(lapply(row, as.numeric), length(score_sql)), groups = as.numeric(row[["n_groups"]]))
}

# Fisher information, score meat and fixest-compatible GLM corrections.
dbglm_ssc = function(ssc) {
  if (is.character(ssc)) {
    if (length(ssc) != 1L || is.na(ssc) || !ssc %in% c("full", "nested", "none")) {
      stop("`ssc` must be full, nested, none, or a correction-policy list.", call. = FALSE)
    }
    ssc = list(K.fixef = switch(ssc,
      full = "full",
      nested = "nonnested",
      none = "none"
    ))
  }
  defaults = list(
    K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE,
    G.adj = TRUE, G.df = "min", t.df = "min"
  )
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
  numeric_df = is.numeric(defaults[["t.df"]])
  if (numeric_df && (length(defaults[["t.df"]]) != 1L || !is.finite(defaults[["t.df"]]) || defaults[["t.df"]] < 1)) {
    stop("Numeric t.df must be at least one.", call. = FALSE)
  }
  for (name in c("K.fixef", "G.df", if (!numeric_df) "t.df")) {
    if (!is.character(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` policies must be single character strings.", call. = FALSE)
    }
  }
  if (defaults[["K.exact"]] || !defaults[["K.fixef"]] %in% c("full", "none", "nonnested") ||
    !defaults[["G.df"]] %in% c("min", "conventional") || !numeric_df && !defaults[["t.df"]] %in% c("min", "conventional")) {
    stop("Unsupported `ssc` settings.", call. = FALSE)
  }
  defaults
}

dbglm_covariance = function(conn, table, x, fit, family, kind, ssc, n, cell = FALSE) {
  inverse = dbglm_solve(fit[["H"]])[["inverse"]]
  p = length(x)
  if (ssc[["K.adj"]] && n <= p) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  if (kind == "iid") {
    V = inverse
    if (ssc[["K.adj"]]) V = V * (n - 1) / (n - p)
    return(list(matrix = V, clusters = NULL))
  }
  cte = dbglm_working_sql(table, x, fit[["beta"]], family[["link"]], cell = cell)
  if (kind == "hc1") {
    row = backend_get_query(conn, paste0(
      cte, "SELECT ",
      paste(dbglm_moment_terms(x, if (cell) "score_sq" else "score*score"), collapse = ", "), " FROM working"
    ))
    G = NULL
  } else {
    scores = sprintf("SUM((%s)*score) AS s_%d", x, seq_along(x))
    terms = dbglm_moment_terms(paste0("s_", seq_along(x)), "1.0")
    row = backend_get_query(conn, paste0(
      cte,
      ", cluster_scores AS (SELECT __dbglm_cluster, ", paste(scores, collapse = ", "),
      " FROM working GROUP BY __dbglm_cluster) SELECT COUNT(*) AS n_clusters, ",
      paste(terms, collapse = ", "), " FROM cluster_scores"
    ))
    G = as.numeric(row[["n_clusters"]])
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  meat = dbglm_read_information(row, p)
  V = inverse %*% meat %*% inverse
  if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - p)
  if (kind == "cluster" && ssc[["G.adj"]]) V = V * G / (G - 1)
  list(matrix = V, clusters = G)
}

dbglm_fe_covariance = function(conn, x, fit, kind, ssc, n, counts, cell = FALSE) {
  p = length(x)
  K = if (kind == "cluster") counts[["covariance"]] else counts[["covariance"]]
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  inverse = if (p) dbglm_solve(fit[["H"]])[["inverse"]] else matrix(numeric(), 0, 0)
  G = NULL
  if (kind == "cluster" && !p) {
    G = as.numeric(backend_get_query(conn, paste("SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM", fit[["working"]]))[["n"]])
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  if (!p) return(list(matrix = inverse, clusters = G, K = K))
  if (kind == "iid") V = inverse else {
    if (!is.null(fit[["raw_projected"]])) {
      from = paste0(" FROM ", fit[["raw_projected"]], " b")
      x_cols = paste0("b.", x, "_tilde")
    } else {
      from = paste0(" FROM ", fit[["working"]], " b JOIN ", fit[["projected"]], " p USING (__dbglm_key)")
      x_cols = paste0("p.", x)
    }
    if (kind == "hc1") {
      sql = paste0("SELECT ", paste(dbglm_moment_terms(x_cols, if (cell) "b.score_sq" else "b.score*b.score"), collapse = ","), from)
      row = backend_get_query(conn, sql)
    } else {
      scores = sprintf("SUM(%s*b.score) AS s%d", x_cols, seq_along(x))
      sql = paste0(
        "SELECT COUNT(*) AS __n_clusters,", paste(dbglm_moment_terms(paste0("s", seq_along(x)), "1.0"), collapse = ","), " FROM (SELECT ",
        paste(scores, collapse = ","), from, " GROUP BY __dbglm_cluster) t"
      )
      row = backend_get_query(conn, sql)
      G = as.numeric(row[["__n_clusters"]])
      if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    }
    meat = dbglm_read_information(row, p)
    V = inverse %*% meat %*% inverse
  }
  if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - K)
  if (kind == "cluster" && ssc[["G.adj"]]) V = V * G / (G - 1)
  list(matrix = V, clusters = G, K = K)
}
