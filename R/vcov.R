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
  sandwich_bread = NULL
) {
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
    scale_hc1 = nobs_orig / df_res
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
    scale_cr1 = (n_clusters / (n_clusters - 1)) * ((nobs_orig - 1) / (nobs_orig - n_params))
    vcov_mat = scale_cr1 * (sandwich_bread %*% meat %*% t(sandwich_bread))
    attr(vcov_mat, "type") = "cluster"
    attr(vcov_mat, "n_clusters") = n_clusters
  } else {
    # IID case (same for all strategies)
    sigma2 = rss / df_res
    vcov_mat = sigma2 * XtX_inv
    attr(vcov_mat, "type") = "iid"
  }
  attr(vcov_mat, "n_params") = if (is.null(n_params)) ncol(XtX_inv) else n_params
  attr(vcov_mat, "df.t") = if (vcov_type == "cluster") min(n_clusters)-1 else df_res
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
    cross_terms = unlist(lapply(seq_along(resid), function(r) sprintf("SUM(score_%d_%d * score_%d_%d) AS meat_%d_%d_%d",
      r, pairs[, 1], r, pairs[, 2], r, pairs[, 1], pairs[, 2])))
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
