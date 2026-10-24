#' Methods for database binomial and Poisson GLMs
#' @param object,x A `dbglm` model.
#' @param ... Additional method arguments are not supported and produce an error.
#'   New-row prediction with absorbed effects requires `retain = "prediction"`
#'   when fitting. Unknown or non-estimable FE combinations return `NA`.
#' @name dbglm-methods
#' @importFrom stats logLik nobs df.residual predict
NULL

#' @rdname dbglm-methods
#' @param type For `logLik`, `"fixest"` reports the likelihood using the
#'   `fixest::feglm()` convention (the default); `"weighted"` returns the
#'   objective with unrounded observation weights. These differ for fractional
#'   binomial weights. For `predict`, `"link"` returns the linear predictor and
#'   `"response"` returns the response mean.
#' @export
logLik.dbglm = function(object, type = c("fixest", "weighted"), ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  value = if (type == "weighted" && !is.null(object[["weighted_loglik"]])) object[["weighted_loglik"]] else object[["loglik"]]
  structure(value, class = "logLik", df = if (isFALSE(object[["finite_mle"]])) NA_integer_ else object[["n_params_model"]],
    nobs = object[["nobs_orig"]], finite_mle = object[["finite_mle"]])
}

#' @rdname dbglm-methods
#' @export
nobs.dbglm = function(object, ...) object[["nobs_orig"]]

#' @rdname dbglm-methods
#' @export
df.residual.dbglm = function(object, ...) object[["df_residual"]]

#' @rdname dbglm-methods
#' @param parm Coefficient names or positions, defaulting to all coefficients.
#' @param level Confidence level.
#' @param fe Must be `FALSE` for GLM confidence intervals. This argument is
#'   accepted for shared tidying; absorbed-effect intervals are not supported.
#' @export
confint.dbglm = function(object, parm, level = 0.95, fe = FALSE, ...) {
  check_unused_args(list(...))
  if (!identical(fe, FALSE)) stop("Absorbed-effect confidence intervals are not supported; use `fe = FALSE`.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must lie strictly between zero and one.", call. = FALSE)
  }
  beta = stats::coef(object)
  if (missing(parm)) parm = names(beta)
  if (is.numeric(parm)) parm = names(beta)[parm]
  if (anyNA(parm) || any(!parm %in% names(beta))) stop("Unknown coefficient in `parm`.", call. = FALSE)
  critical = stats::qnorm(c((1 - level) / 2, (1 + level) / 2))
  result = outer(sqrt(diag(stats::vcov(object)))[parm], critical) + beta[parm]
  colnames(result) = paste0(format(100 * c((1 - level) / 2, (1 + level) / 2)), " %")
  result
}

#' @rdname dbglm-methods
#' @param newdata Data frame containing regressors, fitted FE keys and varying
#'   slope variables, and any fitted offset column. Outcomes and fitting weights
#'   are not required. Missing inputs preserve their positions as `NA` predictions.
#' @export
predict.dbglm = function(object, newdata, type = c("link", "response"), ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  if (missing(newdata)) stop("`newdata` is required; dbglm does not retain input rows.", call. = FALSE)
  if (length(object[["fe"]]) && is.null(object[["fitted_state"]])) {
    stop("New-row prediction for absorbed GLMs requires fitting with retain='prediction'.", call. = FALSE)
  }
  newdata = as.data.frame(newdata)
  for (a in names(object[["factor_sources"]])) {
    source = object[["factor_sources"]][[a]]
    if (!source %in% names(newdata)) stop("newdata is missing columns: ", source, call. = FALSE)
    newdata[[a]] = newdata[[source]]
  }
  state = object[["fitted_state"]]
  required = unique(c(object[["model_spec"]][["source_vars"]], object[["offset"]],
    unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(object[["varying_slopes"]])))
  absent = setdiff(required, names(newdata))
  if (length(absent)) stop("newdata is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  valid = rep(TRUE, nrow(newdata))
  for (name in required) {
    valid = valid & !is.na(newdata[[name]])
    # Nonmissing integer64 values are finite; avoid a lossy conversion of
    # categorical keys merely to test finiteness.
    if (is.numeric(newdata[[name]]) && !inherits(newdata[[name]], "integer64")) {
      valid = valid & is.finite(newdata[[name]])
    }
  }
  prediction = rep(NA_real_, nrow(newdata))
  if (!any(valid)) return(prediction)
  data = newdata[valid, , drop = FALSE]
  positions = which(valid)
  groups = NULL
  if (!is.null(state)) {
    groups = prediction_match_groups(state, data)
    known = if (length(groups)) Reduce(`&`, lapply(groups, function(group) !is.na(group))) else rep(TRUE, nrow(data))
    if (!any(known)) return(prediction)
    if (!all(known)) {
      data = data[known, , drop = FALSE]
      groups = lapply(groups, function(group) group[known])
      positions = positions[known]
    }
  }
  design = as.matrix(evaluate_model_spec(object[["model_spec"]], data))
  colnames(design) = object[["model_spec"]][["columns"]][["name"]]
  if (object[["has_intercept"]]) design = cbind("(Intercept)" = 1, design)
  beta = stats::coef(object)
  names(beta) = object[["design_names"]]
  eta = if (!is.null(state)) predict_from_fitted_state(state, design, data, groups = groups) else
    as.vector(design[, names(beta), drop = FALSE] %*% beta)
  if (!is.null(object[["offset"]])) eta = eta + as.numeric(data[[object[["offset"]]]])
  prediction[positions] = if (type == "link") eta else object[["family"]][["linkinv"]](eta)
  prediction
}

#' @rdname dbglm-methods
#' @param conf.int Include confidence intervals?
#' @param conf.level Confidence level used by `tidy`.
#' @param collin Include structurally collinear coefficients as NA?
#' @export
tidy.dbglm = function(x, conf.int = FALSE, conf.level = 0.95, collin = FALSE, ...) {
  tidy.dbreg(x, conf.int = conf.int, conf.level = conf.level, collin = collin, ...)
}

#' @rdname dbglm-methods
#' @export
glance.dbglm = function(x, ...) {
  data.frame(nobs = x[["nobs_orig"]], nobs.input = x[["nobs_input"]], nobs.dropped = x[["nobs_dropped"]],
    logLik = x[["loglik"]], deviance = x[["deviance"]], AIC = stats::AIC(x), BIC = stats::BIC(x),
    df.residual = x[["df_residual"]], n.params = x[["n_params_model"]], family = x[["family"]][["family"]],
    link = x[["family"]][["link"]], vcov.type = x[["vcov_type"]], converged = x[["converged"]],
    finite.mle = x[["finite_mle"]], iterations = x[["iterations"]], score.max = x[["score_max"]])
}

#' @rdname dbglm-methods
#' @export
print.dbglm = function(x, ...) {
  cat("Database ", x[["family"]][["family"]], " GLM (", x[["family"]][["link"]], "), ", x[["nobs_orig"]],
    " observations; ", x[["iterations"]], " IRLS iterations\n", sep = "")
  if (isFALSE(x[["finite_mle"]])) cat("Boundary probability approximation; no finite MLE established; inference unavailable.\n")
  print(x[["coeftable"]])
  invisible(x)
}
