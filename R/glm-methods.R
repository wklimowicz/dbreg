#' Methods for database binomial and Poisson GLMs
#' @param object,x A `dbglm` model.
#' @param ... Additional method arguments. New-row prediction with absorbed effects is not supported.
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
#' @export
confint.dbglm = function(object, parm, level = 0.95, ...) {
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
#' @param newdata Data frame containing regressors and any fitted offset column.
#' @export
predict.dbglm = function(object, newdata, type = c("link", "response"), ...) {
  type = match.arg(type)
  if (missing(newdata)) stop("`newdata` is required; dbglm does not retain input rows.", call. = FALSE)
  if (length(object[["fe"]])) stop("New-row prediction for absorbed GLMs is not implemented.", call. = FALSE)
  newdata = as.data.frame(newdata)
  for (a in names(object[["factor_sources"]])) {
    source = object[["factor_sources"]][[a]]
    if (!source %in% names(newdata)) stop("newdata is missing columns: ", source, call. = FALSE)
    newdata[[a]] = newdata[[source]]
  }
  required = unique(c(object[["model_spec"]][["source_vars"]], object[["offset"]]))
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
  design = as.matrix(evaluate_model_spec(object[["model_spec"]], data))
  colnames(design) = object[["model_spec"]][["columns"]][["name"]]
  if (object[["has_intercept"]]) design = cbind("(Intercept)" = 1, design)
  beta = stats::coef(object)
  names(beta) = object[["design_names"]]
  eta = as.vector(design[, names(beta), drop = FALSE] %*% beta)
  if (!is.null(object[["offset"]])) eta = eta + as.numeric(data[[object[["offset"]]]])
  prediction[valid] = if (type == "link") eta else object[["family"]][["linkinv"]](eta)
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
