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
#' @param collin Include regressors removed for collinearity with `NA` values?
#' @param ... Additional arguments to tidying method. Currently unused except to
#' handle superseded arguments.
#' @examples
#' mod = dbreg(Temp ~ Wind | Month, data = airquality)
#' tidy(mod, conf.int = TRUE)
#' tidy(mod, conf.int = TRUE, fe = TRUE)
#' glance(mod)
#' @name tidiers
#' @rdname tidiers
#' @export
tidy.dbreg = function(x, conf.int = FALSE, conf.level = 0.95, fe = FALSE,
                       collin = FALSE, ...) {
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
  ct = x[["coeftable"]]
  if (!isTRUE(fe) && !is.null(x[["fe"]])) {
    coef_names = x[["coef_names"]]
    if (!is.null(coef_names)) {
      ct = ct[coef_names, , drop = FALSE]
    } else {
      xvars = x[["xvars"]]
      ct = ct[xvars, , drop = FALSE]
    }
  }
  if (isTRUE(collin)) {
    terms = names(coef(x, fe = fe, collin = TRUE))
    out = data.frame(
      term = terms,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      row.names = NULL
    )
    matched = match(terms, rownames(ct))
    present = !is.na(matched)
    out[present, colnames(ct)] = ct[matched[present], , drop = FALSE]
  } else {
    out = data.frame(term = rownames(ct), ct, row.names = NULL)
  }
  
  if (isTRUE(conf.int)) {
    ci = confint(x, level = conf.level, fe = fe)
    ci_match = match(out[["term"]], rownames(ci))
    out[["conf.low"]] = ifelse(is.na(ci_match), NA_real_, ci[ci_match, 1])
    out[["conf.high"]] = ifelse(is.na(ci_match), NA_real_, ci[ci_match, 2])
  }
  
  return(out)
}

#' glance results
#' @export
generics::glance

#' @importFrom generics glance
#' @rdname tidiers
#' @export
glance.dbreg = function(x, ...) {
  gof_vals = gof(x)
  fe_levels = x[["n_fe_levels"]]
  if (is.null(fe_levels)) {
    fe_levels_text = NA_character_
  } else {
    fe_names = x[["fe_labels"]]
    if (is.null(fe_names)) fe_names = x[["fe"]]
    if (is.null(names(fe_levels)) && length(fe_names) == length(fe_levels)) {
      names(fe_levels) = fe_names
    }
    fe_levels_text = paste(paste0(names(fe_levels), "=", fe_levels), collapse = ";")
  }
  data.frame(
    r.squared = gof_vals["r2"],
    adj.r.squared = gof_vals["adj_r2"],
    rmse = gof_vals["rmse"],
    nobs = x[["nobs_orig"]],
    df.residual = x[["df_residual"]],
    nobs.input = x[["nobs_input"]],
    nobs.dropped = x[["nobs_dropped"]],
    vcov.type = x[["vcov_type"]],
    cluster = if (is.null(x[["cluster_var"]])) NA_character_ else x[["cluster_var"]],
    n.clusters = if (is.null(x[["n_clusters"]])) NA_integer_ else x[["n_clusters"]],
    ssc = if (is.null(x[["ssc"]])) NA_character_ else x[["ssc"]],
    n.params = if (is.null(x[["n_params"]])) NA_integer_ else x[["n_params"]],
    n.fe = length(x[["fe"]]),
    fe.levels = fe_levels_text,
    strategy = x[["strategy"]],
    ap.converged = x[["ap_converged"]],
    ap.iterations = if (is.null(x[["ap_iterations"]])) NA_integer_ else x[["ap_iterations"]],
    ap.tolerance = x[["ap_tolerance"]],
    ap.max.abs = if (is.null(x[["ap_max_abs"]])) NA_real_ else x[["ap_max_abs"]],
    row.names = NULL
  )
}
