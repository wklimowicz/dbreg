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

#' Materialize fixed-effect columns for prediction
#' @keywords internal
materialize_prediction_fes = function(object, newdata) {
  fe = object[["fe"]]
  if (!length(fe)) return(newdata)

  for (i in seq_along(fe)) {
    components = object[["fe_components"]][[i]]
    missing_components = setdiff(components, names(newdata))
    if (length(missing_components)) {
      # Stored compressed data from older objects may only contain the already
      # materialized simple FE column.
      if (length(components) == 1L && fe[[i]] %in% names(newdata)) next
      stop(
        "Fixed effect `", object[["fe_labels"]][[i]],
        "` requires columns: ", paste(missing_components, collapse = ", "),
        call. = FALSE
      )
    }
    values = if (length(components) > 1L) {
      composite_fe_key(newdata, components)
    } else {
      newdata[[components]]
    }
    training_levels = object[["fe_levels"]][[fe[[i]]]]
    if (length(training_levels)) {
      unseen = setdiff(unique(as.character(values)), training_levels)
      unseen = unseen[!is.na(unseen)]
      if (length(unseen)) {
        stop(
          "newdata contains fixed-effect levels not seen during fitting for `",
          object[["fe_labels"]][[i]], "`.",
          call. = FALSE
        )
      }
      newdata[[fe[[i]]]] = factor(values, levels = training_levels)
    } else {
      newdata[[fe[[i]]]] = factor(values)
    }
  }
  newdata
}

weighted_group_mean = function(x, group, weights = NULL) {
  if (is.null(weights)) {
    return(ave(x, group, FUN = function(value) mean(value, na.rm = TRUE)))
  }
  ave(seq_along(x), group, FUN = function(index) {
    sum(x[index] * weights[index], na.rm = TRUE) /
      sum(weights[index], na.rm = TRUE)
  })
}

demean_prediction_data = function(matrix, groups, weights = NULL,
                                   tol = 1e-10, max_iter = 1000L) {
  out = as.matrix(matrix)
  for (iter in seq_len(max_iter)) {
    max_mean = 0
    for (group in groups) {
      means = vapply(seq_len(ncol(out)), function(j) {
        weighted_group_mean(out[, j], group, weights)
      }, numeric(nrow(out)))
      out = out - means
      max_mean = max(max_mean, abs(means), na.rm = TRUE)
    }
    if (is.finite(max_mean) && max_mean <= tol) return(out)
  }
  stop("Prediction demeaning did not converge.", call. = FALSE)
}

dbreg_prediction_matrix = function(object, newdata) {
  spec = object[["model_spec"]]
  if (is.null(spec)) {
    stop("This fitted object does not contain reusable design metadata.", call. = FALSE)
  }
  if (!is.null(object[["weights"]]) && !object[["weights"]] %in% names(newdata)) {
    stop(
      "newdata is missing the fitted weight column: ", object[["weights"]],
      call. = FALSE
    )
  }
  newdata = materialize_prediction_fes(object, newdata)
  design_alias = evaluate_model_spec(spec, newdata)
  design = design_alias
  colnames(design) = spec[["columns"]][["name"]]
  strategy = object[["strategy"]]

  if (strategy == "demean") {
    groups = lapply(object[["fe"]], function(fe) newdata[[fe]])
    weights = if (!is.null(object[["weights"]])) newdata[[object[["weights"]]]] else NULL
    tolerance = object[["ap_tolerance"]]
    if (is.null(tolerance) || is.na(tolerance)) tolerance = 1e-10
    iterations = object[["ap_max_iter"]]
    if (is.null(iterations) || is.na(iterations)) iterations = 1000L
    design = demean_prediction_data(
      design, groups, weights,
      tol = tolerance,
      max_iter = iterations
    )
    y_offset = 0
    if (object[["yvar"]] %in% names(newdata)) {
      y = matrix(newdata[[object[["yvar"]]]], ncol = 1)
      y_within = demean_prediction_data(
        y, groups, weights,
        tol = tolerance,
        max_iter = iterations
      )
      y_offset = y[, 1] - y_within[, 1]
    } else {
      message(
        "Outcome '", object[["yvar"]],
        "' not found in newdata. Returning within-group predictions (demean strategy only)."
      )
    }
    matrix = design
  } else if (strategy == "mundlak") {
    plan = object[["mundlak_plan"]]
    means = matrix(0, nrow(design), nrow(plan))
    weights = if (!is.null(object[["weights"]])) newdata[[object[["weights"]]]] else NULL
    for (i in seq_len(nrow(plan))) {
      source = match(plan[["source_alias"]][[i]], spec[["columns"]][["alias"]])
      means[, i] = weighted_group_mean(
        design[, source], newdata[[plan[["fe"]][[i]]]], weights
      )
    }
    colnames(means) = plan[["name"]]
    matrix = cbind("(Intercept)" = 1, design, means)
    y_offset = 0
  } else if (strategy == "compress") {
    matrix_data = as.data.frame(as.matrix(design_alias))
    for (fe in object[["fe"]]) matrix_data[[fe]] = newdata[[fe]]
    matrix = Matrix::sparse.model.matrix(
      reformulate(c(colnames(design_alias), object[["fe"]])),
      matrix_data
    )
    design_idx = match(spec[["columns"]][["alias"]], colnames(matrix))
    colnames(matrix)[design_idx] = spec[["columns"]][["name"]]
    y_offset = 0
  } else {
    matrix = cbind("(Intercept)" = 1, design)
    y_offset = 0
  }

  coefficients = names(coef(object, fe = TRUE))
  missing = setdiff(coefficients, colnames(matrix))
  if (length(missing)) {
    stop(
      "Prediction matrix is missing fitted columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  list(matrix = matrix[, coefficients, drop = FALSE], offset = y_offset)
}

#' Predict method for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param newdata Data frame for predictions. Required for objects estimated by
#'   strategies that do not retain data.
#' @param interval Type of interval: `"none"`, `"confidence"`, or
#'   `"prediction"`.
#' @param level Confidence level for intervals.
#' @param ... Additional arguments, currently unused.
#'
#' Demean predictions use the fitted expanded design and alternating
#' projections. If the outcome is present, its absorbed component is added back;
#' otherwise within-group predictions are returned. Intervals remain unsupported
#' because uncertainty in absorbed fixed effects is unavailable.
#'
#' @seealso [dbreg()]
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

  betas = coef(object, fe = TRUE)
  prediction = dbreg_prediction_matrix(object, newdata)
  mm = prediction[["matrix"]]

  # Generate predictions
  fit = as.vector(mm %*% betas)
  fit = fit + prediction[["offset"]]
 
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
  metadata = attributes(out)[setdiff(names(attributes(out)), c("dim", "dimnames"))]
  normalized_names = rownames(out)
  dimnames(out) = list(normalized_names, normalized_names)
  if (!isTRUE(fe)) {
    visible = names(coef(object, fe = FALSE))
    visible = intersect(visible, rownames(out))
    out = out[visible, visible, drop = FALSE]
  }
  for (name in names(metadata)) attr(out, name) = metadata[[name]]
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
