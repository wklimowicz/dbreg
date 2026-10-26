# Estimators on sufficient statistics
#
# Strategies aggregate a single cross-moment matrix S in SQL; estimators then
# work purely on (sub-blocks of) S. Each fit_*() returns the same shape:
#   coef, bread, rss, x (regressors kept), scores (score columns for the
#   sandwich meat), project_meat (maps score-space meat to coefficient space)
# so the score pass and vcov scaling are shared by every estimator.
# Centered moments carry location/weight metadata. bread is the inverse
# information in source coordinates; sandwich_bread maps centered score
# coordinates to source coefficients and need not be symmetric. Score passes
# keep their centered columns until the final covariance transformation.

#' Assemble a cross-moment matrix from aggregated SQL sums
#'
#' Reads the `sum_w*` columns produced by [build_weighted_moment_terms()] into
#' one symmetric matrix over `c("(Intercept)"?, vars, y)`. Any (sub-)regression
#' can then be solved by indexing, e.g. `S[x, x]` and `S[x, y]`.
#'
#' @param moment_df One-row data frame of aggregated moments.
#' @param vars Variable names, in the order passed to
#'   [build_weighted_moment_terms()].
#' @param y Name for the outcome row/column.
#' @param intercept Logical. Include the constant (requires the linear sums).
#' @keywords internal
cross_moments = function(moment_df, vars, y, intercept = TRUE) {
  # DBI drivers may return BIGINT aggregates as bit64::integer64. Assignment
  # into a double matrix strips that class without converting its bit pattern.
  # Convert each aggregate before it crosses the SQL-to-linear-algebra boundary.
  moment_df = lapply(moment_df, as.numeric)
  nms = c(if (intercept) "(Intercept)", vars, y)
  S = matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
  S[y, y] = moment_df[["sum_wy_sq"]]
  for (i in seq_along(vars)) {
    S[vars[i], y] = S[y, vars[i]] = moment_df[[sprintf("sum_w%d_y", i)]]
    for (j in seq_len(i)) {
      S[vars[i], vars[j]] = S[vars[j], vars[i]] = moment_df[[sprintf("sum_w%d_%d", j, i)]]
    }
  }
  if (intercept) {
    S[1, 1] = moment_df[["sum_w"]]
    S[1, y] = S[y, 1] = moment_df[["sum_wy"]]
    for (i in seq_along(vars)) {
      S[1, vars[i]] = S[vars[i], 1] = moment_df[[sprintf("sum_w%d", i)]]
    }
  }
  centers = paste0("dbreg_moment_center_", seq_along(c(vars, y)))
  original = paste0("dbreg_moment_original_", seq_along(vars))
  if (length(vars) && all(original %in% names(moment_df))) {
    maximum = vapply(original, function(name) moment_df[[name]], numeric(1))
    floor = moment_df[["sum_w"]] * (128 * .Machine[["double.eps"]] * maximum)^2
    absorbed = vars[diag(S[vars, vars, drop = FALSE]) <= floor]
    S[absorbed, ] = 0
    S[, absorbed] = 0
  }
  if (all(centers %in% names(moment_df))) {
    attr(S, "centers") = setNames(vapply(centers, function(name) moment_df[[name]], numeric(1)), c(vars, y))
    attr(S, "weight_scale") = moment_df[["dbreg_moment_weight_scale"]]
  }
  S
}

# Map centered coefficients back to source coordinates without reconstructing
# the ill-conditioned raw Gram matrix. The same map applies to OLS and IV.
moment_coordinates = function(S, names) {
  A = diag(length(names))
  dimnames(A) = list(names, names)
  centers = attr(S, "centers")
  if ("(Intercept)" %in% names && !is.null(centers)) {
    slopes = setdiff(names, "(Intercept)")
    A["(Intercept)", slopes] = -centers[slopes]
  }
  A
}

restore_moment_fit = function(fit, S, y) {
  A = moment_coordinates(S, fit[["x"]])
  scale = attr(S, "weight_scale")
  if (is.null(scale)) scale = 1
  fit[["moment_coef"]] = fit[["coef"]]
  fit[["centers"]] = attr(S, "centers")
  fit[["sandwich_bread"]] = (A %*% fit[["bread"]]) / scale
  fit[["coef"]] = A %*% fit[["coef"]]
  if ("(Intercept)" %in% fit[["x"]] && !is.null(attr(S, "centers")))
    fit[["coef"]]["(Intercept)", 1L] = fit[["coef"]]["(Intercept)", 1L] + attr(S, "centers")[[y]]
  fit[["bread"]] = (A %*% fit[["bread"]] %*% t(A)) / scale
  fit[["rss"]] = fit[["rss"]] * scale
  fit[["weight_scale"]] = scale
  fit
}

# RSS from moments: y'y - 2b'X'y + b'X'Xb
moment_rss = function(S, y, x, beta) {
  as.numeric(S[y, y] - 2 * crossprod(beta, S[x, y]) + crossprod(beta, S[x, x] %*% beta))
}

# TSS from moments: centered if S carries an intercept, else uncentered
moment_tss = function(S, y) {
  total = S[y, y]
  if ("(Intercept)" %in% rownames(S)) total = total - S["(Intercept)", y]^2 / S["(Intercept)", "(Intercept)"]
  scale = attr(S, "weight_scale")
  if (is.null(scale)) scale = 1
  total * scale
}

#' OLS on sufficient statistics
#'
#' @param S Cross-moment matrix from [cross_moments()].
#' @param y Outcome name.
#' @param x Regressor names.
#' @param drop_collinear Logical. Drop collinear regressors first (via
#'   [detect_collinearity()]).
#' @param verbose Passed to [detect_collinearity()].
#' @keywords internal
fit_ols = function(S, y, x, drop_collinear = FALSE, verbose = FALSE) {
  dropped = character(0)
  if (drop_collinear) {
    collin = detect_collinearity(S[x, x, drop = FALSE], S[x, y, drop = FALSE], verbose = verbose)
    x = collin[["keep_names"]]
    dropped = collin[["drop_names"]]
  }
  sol = solve_with_fallback(S[x, x, drop = FALSE], S[x, y, drop = FALSE])
  beta = sol[["betahat"]]
  rownames(beta) = x
  fit = list(
    coef = beta,
    bread = sol[["XtX_inv"]],
    rss = moment_rss(S, y, x, beta),
    x = x,
    dropped = dropped,
    scores = x,
    project_meat = function(meat) meat
  )
  restore_moment_fit(fit, S, y)
}

#' Linear GMM on sufficient statistics
#'
#' `b = (X'Z W Z'X)^-1 X'Z W Z'y` for a weight matrix `W`. With
#' `W = (Z'Z)^-1` this is 2SLS ([fit_2sls()]); with `W` the inverse of the
#' moment covariance it is (two-step) efficient GMM. The sandwich
#' `A^-1 [P Omega P'] A^-1`, with `P = X'Z W`, has the same form as OLS, with
#' instrument scores projected into coefficient space by `project_meat`.
#'
#' @inheritParams fit_ols
#' @param z Instrument names (included exogenous plus excluded instruments).
#' @param W Weight matrix over `z`.
#' @keywords internal
fit_gmm = function(S, y, x, z, W) {
  S_zx = S[z, x, drop = FALSE]
  P = t(S_zx) %*% W
  A = P %*% S_zx
  if (length(detect_collinearity(A, S[x, y, drop = FALSE])[["drop_names"]])) {
    stop("IV system is not identified or regressors are collinear after instrumentation.")
  }
  sol = solve_with_fallback(A, P %*% S[z, y, drop = FALSE])
  beta = sol[["betahat"]]
  rownames(beta) = x
  fit = list(
    coef = beta,
    bread = sol[["XtX_inv"]],
    rss = moment_rss(S, y, x, beta),
    x = x,
    dropped = character(0),
    scores = z,
    project_meat = function(meat) {
      if (is.null(meat)) {
        return(NULL)
      }
      out = P %*% meat %*% t(P)
      attr(out, "n_clusters") = attr(meat, "n_clusters")
      attr(out, "weight_scale") = attr(meat, "weight_scale")
      out
    }
  )
  restore_moment_fit(fit, S, y)
}

#' 2SLS on sufficient statistics
#'
#' Linear GMM ([fit_gmm()]) with weight matrix `(Z'Z)^-1`.
#'
#' @inheritParams fit_gmm
#' @keywords internal
fit_2sls = function(S, y, x, z) {
  S_zz = S[z, z, drop = FALSE]
  if (length(detect_collinearity(S_zz, S[z, y, drop = FALSE])[["drop_names"]])) {
    stop("Instrument matrix is rank deficient.")
  }
  fit_gmm(S, y, x, z, W = solve_with_fallback(S_zz, diag(ncol(S_zz)))[["XtX_inv"]])
}

#' Score-space meat matrices for fitted models
#'
#' Thin wrapper around [compute_score_meat()] that builds each fit's residual
#' from its coefficients and regressors, and scores all fits in one pass.
#' The fits must share the same score columns (e.g. a 2SLS fit and its first
#' stages, which all use the instruments).
#'
#' @param fits List of [fit_ols()] / [fit_2sls()] results.
#' @param frame List with `cte_sql`, `cte_name`, `cols` (named character
#'   vector mapping variable names to SQL columns in `cte_name`) and
#'   `weights_expr`.
#' @param ys Outcome names (in `frame$cols`), one per fit.
#' @param inputs Inputs environment (for `conn`, `vcov_type_req`,
#'   `cluster_var`).
#' @return List of meat matrices (or `NULL`s for iid), one per fit.
#' @keywords internal
fit_meats = function(fits, frame, ys, inputs) {
  if (!inputs[["vcov_type_req"]] %in% c("hc1", "cluster")) return(vector("list", length(fits)))
  frame[["weight_scale"]] = fits[[1L]][["weight_scale"]]
  scores = fits[[1]][["scores"]]
  if (!all(vapply(fits, function(f) identical(f[["scores"]], scores), logical(1)))) {
    stop("fit_meats() requires fits that share the same score columns.")
  }
  resid_sql = vapply(seq_along(fits), function(k) {
    fit_residual_sql(fits[[k]], frame, ys[[k]], inputs[["conn"]])
  }, character(1))
  score_cols = setdiff(scores, "(Intercept)")
  meats = compute_score_meat(
    inputs = inputs,
    frame = frame,
    resid_sql = resid_sql,
    score_sql = vapply(score_cols, function(name) fit_column_sql(fits[[1L]], frame, name, inputs[["conn"]]), character(1)),
    score_names = score_cols,
    score_intercept = "(Intercept)" %in% scores
  )
  if (is.null(meats)) vector("list", length(fits)) else unname(meats)
}

fit_column_sql = function(fit, frame, name, conn) {
  expression = frame[["cols"]][[name]]
  center = fit[["centers"]][[name]]
  if (is.null(center)) center = 0
  paste0("(", sql_cast_double(paste0("(", expression, ")"), conn), "-(", dbglm_number(center), "))")
}

fit_residual_sql = function(fit, frame, y, conn) {
  b = fit[["moment_coef"]]
  if (is.null(b)) {
    b = fit[["coef"]]
    if ("(Intercept)" %in% fit[["x"]] && !is.null(fit[["centers"]])) {
      slopes = setdiff(fit[["x"]], "(Intercept)")
      b["(Intercept)", 1L] = b["(Intercept)", 1L] - fit[["centers"]][[y]] + sum(fit[["centers"]][slopes] * b[slopes, 1L])
    }
  }
  slopes = setdiff(fit[["x"]], "(Intercept)")
  paste0("(", fit_column_sql(fit, frame, y, conn), "-", sql_linear_predictor(
    vapply(slopes, function(name) fit_column_sql(fit, frame, name, conn), character(1)), b[slopes, 1L],
    if ("(Intercept)" %in% fit[["x"]]) b["(Intercept)", 1L]
  ), ")")
}

#' @rdname fit_meats
#' @param fit,y A single fit and its outcome name.
fit_meat = function(fit, frame, y, inputs) {
  retain_linear_score_frame(fit, frame, y, inputs)
  fit_meats(list(fit), frame, y, inputs)[[1]]
}

center_compressed_design = function(X, center) {
  columns = which(center != 0)
  if (length(columns)) X[, columns] = sweep(as.matrix(X[, columns, drop = FALSE]), 2L, center[columns], "-")
  X
}

compressed_cross_moments = function(X, dat, y, center_columns = seq_len(ncol(X))) {
  center = setNames(rep(0, ncol(X)), colnames(X))
  center_columns = setdiff(center_columns, match("(Intercept)", colnames(X)))
  if ("(Intercept)" %in% names(center) && length(center_columns)) {
    center[center_columns] = as.numeric(Matrix::crossprod(dat[["sum_w"]], X[, center_columns, drop = FALSE]) / sum(dat[["sum_w"]]))
  }
  C = center_compressed_design(X, center)
  H = as.matrix(Matrix::crossprod(C, C * dat[["sum_w"]]))
  xy = as.numeric(Matrix::crossprod(C, dat[["sum_wy"]]))
  S = rbind(cbind(H, xy), c(xy, sum(dat[["sum_wy_sq"]])))
  dimnames(S) = rep(list(c(colnames(X), y)), 2L)
  attr(S, "centers") = c(center[names(center) != "(Intercept)"], setNames(dat[["dbreg_cell_center__"]][[1L]], y))
  attr(S, "weight_scale") = dat[["dbreg_cell_scale__"]][[1L]]
  S
}
