##
## Stats methods for dbreg objects
##

#' Number of observations used in a database regression
#'
#' @param object A \code{dbreg} object.
#' @param ... Additional arguments (currently unused).
#' @return The number of original observations in the estimation sample,
#'   before aggregation or compression by the selected strategy.
#' @importFrom stats nobs
#' @export
nobs.dbreg = function(object, ...) object[["nobs_orig"]]

#' Extract coefficients from dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Should the fixed effects be included? Default is `FALSE`.
#' @param collin Should regressors removed for collinearity be returned as `NA`?
#'   Default is `FALSE`.
#' @param ... Only the superseded `fes` argument is accepted. Other arguments
#'   produce an error.
#' @examples
#' mod = dbreg(Temp ~ Wind | Month, data = airquality)
#' 
#' # coefficients
#' coef(mod)
#' coef(mod, fe = TRUE)  # include fixed effects
#' 
#' # confidence intervals
#' confint(mod)
#' 
#' # variance-covariance matrix
#' vcov(mod)
#' 
#' # predictions
#' head(predict(mod, newdata = airquality))
#' head(predict(mod, newdata = airquality, interval = "confidence"))
#' @export
coef.dbreg = function(object, fe = FALSE, collin = FALSE, ...) {
  
  # superseded args handled through ...
  dots = list(...)
  check_unused_args(dots, allowed = "fes")
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning('The `fes` argument has been superseded by `fe` (without the "s") and will be deprecated in a future `dbreg` release.\n')

    }
  }
  
  ct = object[["coeftable"]]
  out = if (isTRUE(collin) && !is.null(object[["collin.coef"]])) {
    object[["collin.coef"]]
  } else {
    setNames(ct[, "estimate"], rownames(ct))
  }
  
  if (!isTRUE(fe) && !is.null(object[["fe"]])) {
    coef_names = if (isTRUE(collin)) object[["coef_names_all"]] else object[["coef_names"]]
    if (!is.null(coef_names)) {
      out = out[intersect(coef_names, names(out))]
    } else {
      xvars = object[["xvars"]]
      if (isTRUE(collin)) xvars = unique(c(xvars, object[["collin.var"]]))
      out = out[intersect(xvars, names(out))]
    }
  }
  out
}

#' Predict method for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param newdata Data frame for predictions. Required for objects that were
#'   estimated using the `"demean"`, `"mundlak"` and `"moments"` strategies, since `dbreg`
#'   does not retain any data for these estimations.
#' @param interval Type of interval to compute: `"none"` (default), 
#'   `"confidence"`, or `"prediction"`. Note that `"confidence"` intervals
#'   reflect uncertainty in the estimated mean, while `"prediction"` intervals
#'   additionally account for residual variance. See
#'   \code{\link[stats]{predict.lm}} for details.
#' @param level Confidence level for intervals. Default is 0.95.
#' @param ... Additional arguments are not supported and produce an error.
#' 
#' @section Predicting on "demean" strategy objects:
#' 
#' With `retain = "prediction"`, predictions use learned training effects,
#' preserve newdata order and require neither outcomes nor fitting weights.
#' Unknown FE keys and combinations outside the identified joint design span
#' return `NA`. Retained learned-effect state supports `interval = "none"` only.
#' Small dictionaries are portable; oversized one- or two-FE dictionaries can stay
#' in their owning database session. See [close.dbreg()] for resource lifetime.
#'
#' Without prediction retention, `"demean"` models carry two caveats:
#' 
#' 1. Predictions require group means to transform back to the original scale.
#' If `newdata` contains the outcome variable, group means are computed from
#' `newdata` and used to return level predictions. If the outcome is absent,
#' within-group predictions (deviations from group means) are returned instead,
#' with a message. These are conditional predictions for the supplied population,
#' not out-of-sample predictions using saved training fixed effects. They change
#' when the rows or outcomes in `newdata` change. Multiway effects use alternating
#' weighted projections; weighted models require their weight column.
#' 
#' 2. Confidence/prediction intervals are not supported. A demeaned model cannot
#' account for uncertainty in the fixed-effects (since these were absorbed at
#' estimation time), which in turn would yield intervals that are too narrow.
#' Requesting intervals for `"demean"` strategy models will return point
#' predictions with a message. Users should re-estimate with a different
#' strategy if intervals are needed.
#'
#' @section Prediction design:
#' Factor levels and treatment coding are retained from estimation. Unseen levels
#' produce an error; incomplete rows produce `NA` without changing row order.
#' Compressed fits support prediction on their stored cells without raw outcomes
#' or weights. Mundlak means are computed from the supplied population, using
#' fitted-model weights when present; they are not saved training group means.
#'
#' @seealso [dbreg()] for examples.
#' @inherit coef.dbreg examples
#' @importFrom stats model.matrix reformulate ave update as.formula coef vcov qt
#' @importFrom Matrix sparse.model.matrix
#' @export
predict.dbreg = function(
  object,
  newdata = NULL,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  ...
) {
  check_unused_args(list(...))
  interval = match.arg(interval)
  stored_compression = identical(object[["strategy"]], "compress") &&
    (is.null(newdata) || identical(newdata, object[["data"]]))
  if (!is.null(object[["fitted_state"]]) && !stored_compression) return(predict_linear_state(object, newdata, interval))
  strategy = object[["strategy"]] 

  # Demean strategy doesn't support intervals (FE uncertainty not available)
  if (strategy == "demean" && interval != "none") {
    message(
      "Confidence/prediction intervals not supported for demean strategy ",
      "(fixed effect uncertainty unavailable). Returning point predictions."
    )
    interval = "none"
  }

  stored = is.null(newdata) || (strategy == "compress" && identical(newdata, object[["data"]]))
  if (stored) {
    if (strategy == "compress" && !is.null(object[["data"]])) {
      newdata = object[["data"]]
    } else {
      stop("newdata is required for predictions, as dbreg does not retain the original data.")
    }
  }

  # Extract common components from object
  betas = coef(object, fe = TRUE)
  fml = object[["fml"]]
  fe = object[["fe"]]
  yvar = object[["yvar"]]
  intercept = "(Intercept)" %in% names(betas)

  required = unique(c(all.vars(stats::formula(fml, lhs = 0, rhs = 1)),
    unlist(object[["composite_fes"]]), setdiff(fe, names(object[["composite_fes"]])),
    if (strategy %in% c("demean", "mundlak")) object[["weights"]],
    if (strategy == "demean" && yvar %in% names(newdata)) yvar))
  if (stored) required = names(newdata)
  missing = setdiff(required, names(newdata))
  if (length(missing)) stop("newdata is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  n_new = nrow(newdata)
  valid = Reduce(`&`, lapply(newdata[, required, drop = FALSE], function(x) !is.na(x)))
  if (!any(valid)) {
    blank = rep(NA_real_, n_new)
    return(if (interval == "none") blank else data.frame(fit = blank, lwr = blank, upr = blank))
  }
  newdata = newdata[valid, , drop = FALSE]
  newdata = prediction_factors(object, newdata)

  # Ensure FE columns are factors
  for (f in fe) {
    if (!is.factor(newdata[[f]])) newdata[[f]] = factor(newdata[[f]])
  }

  # Rebuild the design from the model's own terms (so factors, interactions and
  # no-intercept models expand exactly as in estimation), then align its
  # columns to the coefficients by name
  terms_rhs = attr(terms(formula(fml, lhs = 0, rhs = 1)), "term.labels")
  # i() terms: build their columns directly (model.matrix() can't expand them)
  i_cols = object[["i_columns"]]
  if (NROW(i_cols)) {
    terms_rhs = terms_rhs[vapply(terms_rhs, function(t) is.null(parse_i_term(t)), logical(1))]
    for (k in seq_len(nrow(i_cols))) {
      col = as.numeric(as.character(newdata[[i_cols[["var"]][k]]]) == i_cols[["level"]][k])
      if (!is.na(i_cols[["by"]][k])) col = col * newdata[[i_cols[["by"]][k]]]
      newdata[[i_cols[["name"]][k]]] = col
    }
    terms_rhs = c(terms_rhs, paste0("`", i_cols[["name"]], "`"))
  }
  if (strategy %in% c("compress", "moments")) {
    # FE enter as dummies (compress)
    terms_rhs = c(terms_rhs, fe)
  }
  if (!length(terms_rhs)) {
    terms_rhs = "1"
  }
  fml_mm = reformulate(terms_rhs, intercept = intercept)
  if (stored && !is.null(object[["compressed_design"]])) {
    spec = object[["compressed_design"]]
    data = newdata
    names(data) = unname(spec[["source_names"]][names(data)])
    mm = fitted_model_matrix(spec[["formula"]], data)
    mm = mm[, names(spec[["columns"]]), drop = FALSE]
    colnames(mm) = unname(spec[["columns"]])
  } else mm = fitted_model_matrix(fml_mm, newdata)
  if (NROW(i_cols)) {
    is_i = colnames(mm) %in% paste0("`", i_cols[["name"]], "`")
    colnames(mm)[is_i] = gsub("`", "", colnames(mm)[is_i], fixed = TRUE)
  }
  if (strategy == "mundlak") {
    # Group means of each design column, one set per FE
    mm_means = mapply(
      function(term, f) weighted_group_mean(as.numeric(mm[, term]), newdata[[f]],
        if (is.null(object[["weights"]])) NULL else newdata[[object[["weights"]]]]),
      object[["mundlak_means"]][["term"]],
      object[["mundlak_means"]][["fe"]]
    )
    mm_means = matrix(mm_means, nrow = nrow(mm), dimnames = list(NULL, object[["mundlak_means"]][["name"]]))
    mm = cbind(mm, mm_means)
  }
  missing_cols = setdiff(names(betas), colnames(mm))
  if (length(missing_cols)) {
    stop(
      "predict() could not rebuild the design for coefficient(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  mm = mm[, names(betas), drop = FALSE]

  if (strategy == "demean") {
    # demean: sweep the FE out of the regressors (and recover the outcome's
    # group means, if available, to return level predictions)
    has_y = yvar %in% names(newdata)
    groups = lapply(fe, function(f) newdata[[f]])
    weights = if (is.null(object[["weights"]])) NULL else newdata[[object[["weights"]]]]
    sweep_fe = function(v) demean_prediction_data(v, groups, weights)
    mm = sweep_fe(mm)
    if (has_y) {
      y = as.numeric(newdata[[yvar]])
      y_group_mean = y - drop(sweep_fe(matrix(y, ncol = 1L)))
    } else {
      y_group_mean = 0
      message(sprintf(
        "Outcome '%s' not found in newdata. Returning within-group predictions (demean strategy only).",
        yvar
      ))
    }
  }

  # Generate predictions
  fit = as.vector(mm %*% betas)
  if (strategy == "demean") {
    fit = fit + y_group_mean
  }
 
  if (interval != "none") {
    vcovm = vcov(object, fe = TRUE)
    dof = object[["df_residual"]]
    if (interval == "confidence") {
      ses = sqrt(Matrix::rowSums((mm %*% vcovm) * mm))
    } else if (interval == "prediction") {
      # Use stored RSS from training to estimate sigma^2
      rss = attr(object[["vcov"]], "rss")
      if (is.null(rss)) {
        stop("Prediction intervals require RSS, which is not stored in this model.")
      }
      sig2 = rss / dof
      ses = sqrt(Matrix::rowSums((mm %*% vcovm) * mm) + sig2)
    }
    a = (1 - level) / 2
    t_crit = qt(1 - a, dof)
    lwr = fit - t_crit * ses
    upr = fit + t_crit * ses
    fit = data.frame(fit = fit, lwr = lwr, upr = upr)
  }

  if (interval == "none") {
    out = rep(NA_real_, n_new)
    out[valid] = fit
  } else {
    out = data.frame(fit = rep(NA_real_, n_new), lwr = NA_real_, upr = NA_real_)
    out[valid, ] = fit
  }
  return(out)
}

#' Variance-covariance matrix for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Include covariance rows and columns for explicit fixed effects?
#'   Default is `FALSE`, matching [coef()].
#' @param ... Additional covariance requests are not yet supported and produce
#'   an error. This method returns the covariance calculated when fitting.
#' @inherit coef.dbreg examples
#' @export
vcov.dbreg = function(object, fe = FALSE, ...) {
  check_unused_args(list(...))
  out = object[["vcov"]]
  metadata = lapply(c("type", "n_clusters", "n_params", "rss", "tss"), function(nm) attr(out, nm))
  names(metadata) = c("type", "n_clusters", "n_params", "rss", "tss")
  # The table already has public names (including interaction separators).
  # Its rows and the stored covariance share the fitted parameter order.
  dimnames(out) = rep(list(rownames(object[["coeftable"]])), 2L)
  if (!isTRUE(fe) && length(object[["fe"]])) {
    visible = names(coef(object, fe = fe))
    out = out[visible, visible, drop = FALSE]
  }
  for (nm in names(metadata)) if (!is.null(metadata[[nm]])) attr(out, nm) = metadata[[nm]]
  out
}

#' Confidence intervals for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level the confidence level required. Default is 0.95.
#' @param fe Should the fixed effects be included? Default is `FALSE`.
#' @param ... Only the superseded `fes` argument is accepted. Other arguments
#'   produce an error.
#' @inherit coef.dbreg examples
#' @importFrom stats qt
#' @export
confint.dbreg = function(object, parm, level = 0.95, fe = FALSE, ...) {

  # superseded args handled through ...
  dots = list(...)
  check_unused_args(dots, allowed = "fes")
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning('The `fes` argument has been superseded by `fe` (without the "s") and will be deprecated in a future `dbreg` release.\n')

    }
  }

  ct = object[["coeftable"]]
  
  if (!isTRUE(fe) && !is.null(object[["fe"]])) {
    coef_names = object[["coef_names"]]
    if (!is.null(coef_names)) {
      ct = ct[coef_names, , drop = FALSE]
    } else {
      xvars = object[["xvars"]]
      ct = ct[xvars, , drop = FALSE]
    }
  }
  
  cf = ct[, "estimate"]
  ses = ct[, "std.error"]
  df = object[["df_residual"]]
  
  a = (1 - level) / 2
  t_crit = qt(1 - a, df)
  ci = cbind(cf - t_crit * ses, cf + t_crit * ses)
  rownames(ci) = rownames(ct)
  colnames(ci) = sprintf("%.1f %%", 100 * c(a, 1 - a))
  
  if (!missing(parm)) {
    ci = ci[parm, , drop = FALSE]
  }
  
  ci
}

