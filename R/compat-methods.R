#' Predictions from compatibility entry points
#' @param object A model from [feols()] or [feglm()].
#' @param newdata Optional R prediction data.
#' @param type Response (default) or link; identical for linear models.
#' @param offset Optional newdata offset column/formula or R numeric vector.
#'   Required for newdata when fitting used a positional offset vector.
#' @param interval,level Native OLS interval controls.
#' @param conn,table,path,row_id,prediction_table Native SQL prediction arguments.
#' @param sample Estimation rows or original rows with missing omissions.
#' @param ... Unsupported arguments error.
#' @return Native numeric vectors or keyed database predictions.
#' @export
predict.dbreg_fixest = function(object, newdata = NULL, type = c("response", "link"), offset = NULL,
  interval = c("none", "confidence", "prediction"), level = .95,
  conn = NULL, table = NULL, path = NULL, row_id = NULL, prediction_table = NULL,
  sample = c("estimation", "original"), ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  interval = match.arg(interval)
  sample = match.arg(sample)
  sql = any(!vapply(list(conn,table,path,row_id,prediction_table), is.null, logical(1)))
  supplied = !is.null(newdata) || sql
  compatibility = object[["compatibility"]]
  if (isTRUE(compatibility[["lean"]]) && (!supplied || length(object[["fe"]]))) stop("This prediction requires state released by lean=TRUE.", call. = FALSE)
  if (!is.null(offset) && !supplied) stop("A prediction offset requires newdata or a database source.", call. = FALSE)
  if (supplied && isTRUE(compatibility[["vector_offset"]]) && is.null(offset)) stop("Newdata prediction after a vector offset requires an explicit offset.", call. = FALSE)
  if (!is.null(offset)) {
    if (inherits(offset, "formula")) {
      if (length(offset) != 2L || !is.symbol(offset[[2L]])) stop("Prediction offset formulas must name one column.", call. = FALSE)
      offset = as.character(offset[[2L]])
    }
    if (is.numeric(offset)) {
      if (is.null(newdata) || !is.data.frame(newdata) || !(length(offset) %in% c(1L,nrow(newdata)))) stop("Numeric prediction offsets require R newdata of the corresponding size.", call. = FALSE)
      column = safe_identifiers("dbreg_prediction_offset__", taken = names(newdata), fold_case = TRUE)
      newdata[[column]] = if (length(offset) == 1L) rep(offset,nrow(newdata)) else offset
      offset = column
    }
    validate_linear_offset(offset)
    object[["offset"]] = offset
  }
  if (inherits(object, "dbglm")) {
    if (interval != "none") stop("Intervals are unavailable for the shared GLM/varying-slope prediction path.", call. = FALSE)
    return(predict.dbglm(object,newdata,type,conn,table,path,row_id,prediction_table,sample))
  }
  predict.dbreg(object,newdata,interval,level,conn,table,path,row_id,prediction_table,sample)
}

#' Cached compatibility likelihood
#' @param object A compatibility model.
#' @param type GLM likelihood convention, as in [logLik.dbglm()].
#' @param ... Unsupported arguments error.
#' @return A logLik object, also available after lean/release.
#' @export
logLik.dbreg_fixest = function(object, type = "fixest", ...) {
  check_unused_args(list(...))
  if (object[["compatibility"]][["estimator"]] == "feols") {
    if (!identical(type,"fixest")) stop("OLS supports only its standard likelihood convention.", call. = FALSE)
    return(object[["compatibility"]][["loglik"]])
  }
  logLik.dbglm(object,type)
}

#' Residuals from compatibility models
#' @param object A compatibility model with live retained state.
#' @param type Native residual type; working residuals require feglm.
#' @param na.rm Omit rows excluded from estimation?
#' @param ... Unsupported arguments error.
#' @return Native residual vector or keyed relation.
#' @export
residuals.dbreg_fixest = function(object, type = c("response","deviance","pearson","working"), na.rm = TRUE, ...) {
  type = match.arg(type)
  if (type == "working" && object[["compatibility"]][["estimator"]] == "feols") stop("Working residuals require feglm().", call. = FALSE)
  residuals.dbreg(object,type,na.rm,...)
}

#' Information criteria for compatibility models
#' Uses the model parameter count, excluding the additional OLS likelihood df,
#' matching fixest's AIC/BIC reporting convention.
#' @param object A compatibility model.
#' @param ... Additional fitted models; returns a numeric vector like fixest.
#' @param k AIC penalty for the first model.
#' @return Numeric information criterion values.
#' @importFrom stats AIC BIC
#' @export
AIC.dbreg_fixest = function(object, ..., k = 2) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) stop("k must be one finite number.", call. = FALSE)
  others = list(...)
  if (length(others) && any(vapply(others, stats::nobs, numeric(1)) != nobs(object))) warning("Models are not all fitted to the same number of observations.", call. = FALSE)
  c(-2*as.numeric(stats::logLik(object))+k*object[["n_params_model"]], vapply(others, stats::AIC, numeric(1)))
}

#' @rdname AIC.dbreg_fixest
#' @export
BIC.dbreg_fixest = function(object, ...) {
  others = list(...)
  if (length(others) && any(vapply(others, stats::nobs, numeric(1)) != nobs(object))) warning("Models are not all fitted to the same number of observations.", call. = FALSE)
  c(-2*as.numeric(stats::logLik(object))+log(nobs(object))*object[["n_params_model"]], vapply(others, stats::BIC, numeric(1)))
}
