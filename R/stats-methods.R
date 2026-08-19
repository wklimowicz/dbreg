##
## Stats methods for dbreg objects
##

#' Extract coefficients from dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Should the fixed effects be included? Default is `FALSE`.
#' @param collin Should regressors removed for collinearity be returned as `NA`?
#'   Default is `FALSE`.
#' @param ... Additional arguments. Currently unused, except to capture
#' superseded arguments.
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
#'   estimated using the `"mundlak"` and `"moments"` strategies, since `dbreg`
#'   does not retain any data for these estimations.
#' @param interval Type of interval to compute: `"none"` (default), 
#'   `"confidence"`, or `"prediction"`. Note that `"confidence"` intervals
#'   reflect uncertainty in the estimated mean, while `"prediction"` intervals
#'   additionally account for residual variance. See
#'   \code{\link[stats]{predict.lm}} for details.
#' @param level Confidence level for intervals. Default is 0.95.
#' @param ... Additional arguments (currently unused).
#' 
#' @section Predicting on "demean" strategy objects:
#' 
#' Predicting on `dbreg` objects should generally work as expected. However,
#' predictions from `"demean"` strategy models carry two important caveats:
#' 
#' 1. Predictions require group means to transform back to the original scale.
#' If `newdata` contains the outcome variable, group means are computed from
#' `newdata` and used to return level predictions. If the outcome is absent,
#' within-group predictions (deviations from group means) are returned instead,
#' with a message.
#' 
#' 2. Confidence/prediction intervals are not supported. A demeaned model cannot
#' account for uncertainty in the fixed-effects (since these were absorbed at
#' estimation time), which in turn would yield intervals that are too narrow.
#' Requesting intervals for `"demean"` strategy models will return point
#' predictions with a message. Users should re-estimate with a different
#' strategy if intervals are needed.
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
  interval = match.arg(interval)
  strategy = object[["strategy"]] 

  # Demean strategy doesn't support intervals (FE uncertainty not available)
  if (strategy == "demean" && interval != "none") {
    message(
      "Confidence/prediction intervals not supported for demean strategy ",
      "(fixed effect uncertainty unavailable). Returning point predictions."
    )
    interval = "none"
  }

  if (is.null(newdata)) {
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
  xvars = object[["xvars"]]
  yvar = object[["yvar"]]

  factor_refs = object[["factor_refs"]]
  if (length(factor_refs)) {
    for (variable in names(factor_refs)) {
      newdata[[variable]] = stats::relevel(
        factor(newdata[[variable]]),
        ref = as.character(factor_refs[[variable]])
      )
    }
  }

  fe_components = object[["fe_components"]]
  if (!is.null(fe_components)) {
    for (i in seq_along(fe_components)) {
      if (length(fe_components[[i]]) > 1) {
        fe_name = fe[[i]]
        missing_components = setdiff(fe_components[[i]], names(newdata))
        if (length(missing_components)) {
          stop(
            "Composite fixed effect `", object[["fe_labels"]][[i]],
            "` requires columns: ", paste(missing_components, collapse = ", "),
            call. = FALSE
          )
        }
        newdata[[fe_name]] = do.call(
          interaction,
          c(
            as.list(newdata[, fe_components[[i]], drop = FALSE]),
            list(drop = TRUE, lex.order = TRUE)
          )
        )
      }
    }
  }

  # Ensure FE columns are factors
  for (f in fe) {
    newdata[[f]] = factor(newdata[[f]])
  }

  if (strategy == "demean") {
    # demean: compute group means from newdata to demean predictors
    has_y = yvar %in% names(newdata)
    mean_fn = function(y) mean(y, na.rm = TRUE)
    
    if (length(fe) == 1) {
      fe1 = fe[1]
      # Demean X using group means from newdata
      mm = sapply(xvars, function(v) {
        x_mean = ave(newdata[[v]], newdata[[fe1]], FUN = mean_fn)
        newdata[[v]] - x_mean
      })
      
      if (has_y) {
        y_group_mean = ave(newdata[[yvar]], newdata[[fe1]], FUN = mean_fn)
      } else {
        y_group_mean = 0
        message(sprintf(
          "Outcome '%s' not found in newdata. Returning within-group predictions (demean strategy only).",
          yvar
        ))
      }
      
    } else {
      # 2-FE: double demeaning
      fe1 = fe[1]
      fe2 = fe[2]
      
      # Demean X using double-demeaning from newdata
      mm = sapply(xvars, function(v) {
        x_u = ave(newdata[[v]], newdata[[fe1]], FUN = mean_fn)
        x_t = ave(newdata[[v]], newdata[[fe2]], FUN = mean_fn)
        x_o = mean_fn(newdata[[v]])
        newdata[[v]] - x_u - x_t + x_o
      })
      
      if (has_y) {
        y_u = ave(newdata[[yvar]], newdata[[fe1]], FUN = mean_fn)
        y_t = ave(newdata[[yvar]], newdata[[fe2]], FUN = mean_fn)
        y_o = mean_fn(newdata[[yvar]])
        y_group_mean = y_u + y_t - y_o
      } else {
        y_group_mean = 0
        message(sprintf(
          "Outcome '%s' not found in newdata. Returning within-group predictions (demean strategy only).",
          yvar
        ))
      }
    }
    
  } else if (strategy == "mundlak") {
    # Create group means for Mundlak prediction
    gmeans = c()
    for (x in xvars) {
      for (f in fe) {
        demean_x = paste0(x, "_bar_", f)
        gmeans = c(gmeans, demean_x)
        newdata[[demean_x]] = ave(
          newdata[[x]],
          newdata[[f]],
          FUN = function(y) mean(y, na.rm = TRUE)
        )
      }
    }

    # recast as Mundlak formula
    fml_xvars = formula(fml, lhs = 0, rhs = 1)
    fml = update(fml_xvars, as.formula(paste("~ . +", paste(gmeans, collapse = " + "))))

    mm = sparse.model.matrix(fml, data = newdata)
  } else {
    # compress/moments: use sparse model matrix with FE dummies
    mm = sparse.model.matrix(reformulate(c(xvars, fe)), data = newdata)
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

  return(fit)
}

#' Variance-covariance matrix for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Should covariance rows and columns for explicit fixed-effect
#'   coefficients be included? Default is `FALSE`.
#' @param ... Additional arguments (currently unused).
#' @inherit coef.dbreg examples
#' @export
vcov.dbreg = function(object, fe = FALSE, ...) {
  out = object[["vcov"]]
  normalized_names = standardize_coef_names(rownames(out))
  dimnames(out) = list(normalized_names, normalized_names)
  if (!isTRUE(fe)) {
    visible = names(coef(object, fe = FALSE))
    visible = intersect(visible, rownames(out))
    out = out[visible, visible, drop = FALSE]
  }
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
#' @param ... Additional arguments. Currently unused, except to capture
#' superseded arguments.
#' @inherit coef.dbreg examples
#' @importFrom stats qt
#' @export
confint.dbreg = function(object, parm, level = 0.95, fe = FALSE, ...) {

  # superseded args handled through ...
  dots = list(...)
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
