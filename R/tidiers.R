#' tidy results
#' @export
generics::tidy

#' Tidiers for `dbreg` objects
#' @description Provides \code{broom::tidy} and \code{broom::glance} methods for
#' "dbreg" objects. 
#' @importFrom generics tidy
#' @importFrom stats confint
#' @param x a model of class `dbreg` produced by the \code{\link[dbreg]{dbreg}}
#' function.
#' @param conf.int Logical indicating whether to include confidence intervals.
#' Default is `FALSE`.
#' @param conf.level Confidence level for intervals. Default is 0.95.
#' @param fe Should the fixed effects be tidied too? Default is `FALSE`.
#' @param cluster Optional cluster formula or column names for covariance.
#' @param collin Include columns dropped for collinearity as missing rows.
#' @param vcov Optional covariance request or a supplied numeric covariance
#'   matrix. Named axes are aligned to uniquely named coefficients. Duplicate
#'   names require both axes in fitted order.
#' @param ... Additional arguments to tidying method. Currently unused except to
#' handle superseded arguments.
#' @inheritParams summary.dbreg
#' @examples
#' mod = dbreg(Temp ~ Wind | Month, data = airquality)
#' tidy(mod, conf.int = TRUE)
#' tidy(mod, conf.int = TRUE, fe = TRUE)
#' glance(mod)
#' @name tidiers
#' @rdname tidiers
#' @export
tidy.dbreg = function(x, conf.int = FALSE, conf.level = 0.95, fe = FALSE, ...,
                      vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE, collin = FALSE) {
  # superseded args handled through ...
  dots = list(...)
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning(
        'The `fes` argument has been superseded by `fe` (without the "s") and will be deprecated in a future `dbreg` release.\n'
      )
    }
  }
  custom = is.matrix(vcov) || inherits(vcov, "Matrix")
  if (!custom) x = summary.dbreg(x, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
  table = visible_coefficient_table(x, fe)
  df = linear_test_df(x)
  normal = normal_inference(x)
  if (custom) {
    if (!is.null(cluster) || !is.null(ssc) || isTRUE(fe)) {
      stop("A supplied covariance matrix cannot be combined with cluster, ssc, or fe=TRUE.", call. = FALSE)
    }
    supplied_df = attr(vcov, "df.t")
    if (!is.null(supplied_df)) {
      if (!is.numeric(supplied_df) || length(supplied_df) != 1L || is.na(supplied_df) || supplied_df <= 0) {
        stop("Covariance df.t must be one positive number.", call. = FALSE)
      }
      df = supplied_df
    }
    covariance = tidy_covariance_matrix(vcov, rownames(table))
    table = inference_table(table[, "estimate"], covariance, df, normal)
  }
  positions = seq_len(nrow(table))
  labels = rownames(table)
  if (isTRUE(collin) && !is.null(x[["linear_design"]])) {
    fitted = coefficient_view(x, fe)
    full = coefficient_view(x, fe, collin = TRUE)
    positions = match(full[["fitted"]], fitted[["fitted"]])
    labels = full[["labels"]]
  }
  out = data.frame(term = as.character(labels), table[positions, , drop = FALSE], row.names = NULL)
  if (isTRUE(conf.int)) {
    ci = coefficient_limits(table, conf.level, df, normal)
    out[["conf.low"]] = unname(ci[positions, 1L])
    out[["conf.high"]] = unname(ci[positions, 2L])
  }
  out
}

tidy_covariance_matrix = function(V, labels) {
  V = as.matrix(V)
  if (!is.numeric(V) || !identical(dim(V), rep.int(length(labels), 2L)) || any(!is.finite(V))) {
    stop("vcov must be a finite symmetric covariance matrix for the fitted coefficients.", call. = FALSE)
  }
  if (is.null(dimnames(V))) dimnames(V) = list(labels, labels)
  if (!identical(sort(rownames(V)), sort(labels)) || !identical(sort(colnames(V)), sort(labels))) {
    stop("Covariance names must match fitted coefficients.", call. = FALSE)
  }
  if (!identical(rownames(V), labels) || !identical(colnames(V), labels)) {
    if (anyDuplicated(labels)) stop("Duplicate coefficient names require covariance axes in fitted order.", call. = FALSE)
    V = V[labels, labels, drop = FALSE]
  }
  if (!isSymmetric(unname(V)) || any(diag(V) < 0)) {
    stop("vcov must be a finite symmetric covariance matrix for the fitted coefficients.", call. = FALSE)
  }
  V
}

#' glance results
#' @export
generics::glance

#' @importFrom generics glance
#' @rdname tidiers
#' @export
glance.dbreg = function(x, ...) {
  if (inherits(x, "dbglm")) return(glance_glm(x))
  gof_vals = gof(x)
  data.frame(
    r.squared = gof_vals["r2"],
    adj.r.squared = gof_vals["adj_r2"],
    rmse = gof_vals["rmse"],
    nobs = x[["nobs_orig"]],
    df.residual = x[["df_residual"]],
    nobs.input = glance_value(x[["nobs_input"]], NA_real_),
    nobs.dropped = glance_value(x[["nobs_dropped"]], NA_real_),
    vcov.type = glance_value(attr(x[["vcov"]], "type")),
    n.clusters = glance_value(attr(x[["vcov"]], "n_clusters"), NA_integer_),
    n.params = glance_value(attr(x[["vcov"]], "n_params"), NA_integer_),
    n.fe = length(x[["fe"]]), fe.levels = glance_value(x[["n_fe_levels"]]),
    strategy = glance_value(x[["strategy"]]),
    ap.converged = glance_value(x[["ap_converged"]], NA),
    ap.iterations = glance_value(x[["ap_iterations"]], NA_integer_),
    ap.tolerance = glance_value(x[["ap_tolerance"]], NA_real_),
    ap.max.abs = glance_value(x[["ap_max_abs"]], NA_real_),
    row.names = NULL
  )
}


glance_value = function(value, missing = NA_character_) {
  if (!length(value)) return(missing)
  if (!is.null(names(value))) return(paste(paste0(names(value), "=", value), collapse = ";"))
  if (length(value) == 1L) value else paste(value, collapse = ";")
}

glance_glm = function(x) {
  result = data.frame(nobs = nobs(x), nobs.input = x[["nobs_input"]], nobs.dropped = x[["nobs_dropped"]],
    logLik = x[["loglik"]], deviance = x[["deviance"]], AIC = stats::AIC(x), BIC = stats::BIC(x),
    df.residual = x[["df_residual"]], n.params = x[["n_params_model"]], family = x[["family"]][["family"]],
    link = x[["family"]][["link"]], vcov.type = attr(x[["vcov"]], "type"),
    converged = x[["converged"]], finite.mle = x[["finite_mle"]], iterations = x[["iterations"]],
    score.max = x[["score_max"]])
  null = x[["null_diagnostics"]]
  if (!is.null(null)) {
    metrics = dbglm_null_metrics(null)
    for (field in names(metrics)) result[[field]] = metrics[[field]]
    result[["null.reason"]] = glance_value(null[["fit"]][["reason"]])
    result[["null.logLik.reason"]] = glance_value(null[["reported"]][["reason"]])
    result[["weighted.null.logLik.reason"]] = glance_value(null[["weighted"]][["reason"]])
    result[["pseudo.r.squared.reason"]] = glance_value(null[["reported"]][["ratio_reason"]])
    result[["weighted.pseudo.r.squared.reason"]] = glance_value(null[["weighted"]][["ratio_reason"]])
  }
  result
}

#' Covariance for insight and modelsummary
#' @param x A fitted database regression.
#' @param vcov A covariance type or cluster formula.
#' @param vcov_args A named list of cluster, ssc or vcov_fix options.
#' @param ... Unused display arguments.
#' @return The fitted or requested covariance matrix.
#' @exportS3Method insight::get_varcov
get_varcov.dbreg = function(x, vcov = NULL, vcov_args = NULL, ...) {
  if (is.null(vcov_args)) vcov_args = list()
  if (!is.list(vcov_args) || (length(vcov_args) &&
      (is.null(names(vcov_args)) || anyDuplicated(names(vcov_args)) ||
       !all(names(vcov_args) %in% c("cluster", "ssc", "vcov_fix"))))) {
    stop("vcov_args must name cluster, ssc or vcov_fix options.", call. = FALSE)
  }
  if (is.character(vcov) && length(vcov) == 1L && !is.na(vcov)) {
    if (tolower(vcov) == "vcovcl") vcov = if (is.null(vcov_args[["cluster"]])) "cluster" else NULL
    if (identical(tolower(vcov), "vcov")) vcov = NULL
  }
  do.call(vcov.dbreg, c(list(object = x, vcov = vcov), vcov_args))
}
