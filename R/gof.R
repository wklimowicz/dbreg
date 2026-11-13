#' Calculate goodness-of-fit metrics for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param ... Additional arguments (currently unused)
#' @return Named OLS fit statistics or GLM likelihood and deviance statistics.
#' @examples
#' mod = dbreg(Temp ~ Wind | Month, data = airquality)
#' gof(mod)
#' @export
gof = function(object, ...) {
  if (inherits(object, "dbglm")) {
    result = c(logLik = object[["loglik"]], deviance = object[["deviance"]],
      AIC = stats::AIC(object), BIC = stats::BIC(object))
    null = object[["null_diagnostics"]]
    if (!is.null(null)) {
      result = c(result, dbglm_null_metrics(null))
      attr(result, "null_diagnostics") = null
    }
    return(result)
  }
  nobs = object[["nobs_orig"]]
  df_res = object[["df_residual"]]
  
  rss = attr(object[["vcov"]], "rss")
  tss = attr(object[["vcov"]], "tss")
  
  if (is.null(rss) || is.null(tss)) {
    warning("RSS or TSS not available in model object.")
    return(c(r2 = NA_real_, adj_r2 = NA_real_, rmse = NA_real_))
  }
  
  r2_val = 1 - rss / tss
  adj_r2_val = 1 - (rss / df_res) / (tss / (nobs - 1))
  rmse_val = sqrt(rss / nobs)
  
  c(r2 = r2_val, adj_r2 = adj_r2_val, rmse = rmse_val)
}


dbglm_null_metrics = function(null) {
  reported = null[["reported"]]
  weighted = null[["weighted"]]
  value = function(record, field, status) if (isTRUE(record[[status]])) record[[field]] else NA_real_
  c(null.logLik = value(reported, "logLik", "available"),
    weighted.null.logLik = value(weighted, "logLik", "available"),
    pseudo.r.squared = value(reported, "ratio", "ratio_available"),
    weighted.pseudo.r.squared = value(weighted, "ratio", "ratio_available"),
    null.available = as.numeric(null[["fit"]][["available"]]),
    null.logLik.available = as.numeric(reported[["available"]]),
    weighted.null.logLik.available = as.numeric(weighted[["available"]]),
    pseudo.r.squared.available = as.numeric(reported[["ratio_available"]]),
    weighted.pseudo.r.squared.available = as.numeric(weighted[["ratio_available"]]))
}
