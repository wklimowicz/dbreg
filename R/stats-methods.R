#' Number of observations used in a dbreg fit
#'
#' @param object A `dbreg` or `dbivreg` object.
#' @param ... Unused.
#' @return Number of source rows in the estimation sample, not the number of
#'   compressed cells or aggregated moment rows. Weights do not change this count.
#' @importFrom stats nobs
#' @export
nobs.dbreg = function(object, ...) as.numeric(object[["nobs_orig"]])

##
## Stats methods for dbreg objects
##

#' Extract coefficients from dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Should the fixed effects be included? Default is `FALSE`.
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
coef.dbreg = function(object, fe = FALSE, ...) {
  
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
    if (!is.null(object[["linear_design"]])) {
      ct = ct[visible_linear_positions(object, fe), , drop = FALSE]
    } else {
      coef_names = object[["coef_names"]]
      if (!is.null(coef_names)) {
        ct = ct[coef_names, , drop = FALSE]
      } else {
        ct = ct[object[["xvars"]], , drop = FALSE]
      }
    }
  }
  out = ct[, "estimate"]
  names(out) = rownames(ct)
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
#' @param level One finite numeric confidence level strictly between zero and
#'   one. Default is 0.95.
#' @param conn,table,path Optional connection and database prediction source.
#'   Provide exactly one of newdata, table or path.
#' @param row_id Unique nonmissing source key columns for database output.
#' @param prediction_table Optional new session table owned by the caller.
#' @param sample Without `newdata`, a retained fit returns the estimation
#'   sample by default. Use `"original"` to include omitted observations as
#'   missing predictions.
#' @param ... Additional arguments (currently unused).
#'
#' Fitted column construction and category coding are retained independently of
#' public labels. Missing values and unknown fitted categories return `NA` in
#' their input positions. Numeric predictors require numeric or logical input.
#' The original row order is preserved; unrelated columns do not affect it.
#'
#' @section Keyed database output:
#' Supply `row_id` to keep predictions in a database session. This path
#' supports OLS and GLM models on DuckDB and native Databricks. Fixed-effect
#' fits require portable learned effects saved by prediction retention.
#' Mundlak and IV database output are unavailable.
#' Use \code{\link{dbreg-prediction}} to preview, collect or close results.
#'
#' Learned character categories use R's `as.character()` labels for R prediction
#' inputs, preserving the original row keys. Database sources for those categories
#' accept character, exact integer, or logical columns (logical labels are
#' `"TRUE"` and `"FALSE"`); floating-point columns require explicit conversion to
#' the intended character labels. Date categories accept Date columns or canonical
#' ISO date strings; malformed strings produce missing predictions.
#'
#' @section Learned fixed effects:
#' Fits with `retain = "prediction"` or `"all"` use the saved training
#' effects on new R data. Predictions are independent of the other rows
#' supplied and remain available after closing the fit. Unseen groups and
#' non-estimable combinations return missing values. Learned effects support
#' point predictions with `interval = "none"`.
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
  ...,
  conn = NULL,
  table = NULL,
  path = NULL,
  row_id = NULL,
  prediction_table = NULL,
  sample = c("estimation", "original")
) {
  interval = match.arg(interval)
  sample = match.arg(sample)
  sql = sql_prediction_requested(conn, table, path, row_id, prediction_table)
  if (sample == "original" && (!is.null(newdata) || sql)) {
    stop("Original-sample prediction uses the retained fitting source.", call. = FALSE)
  }
  if (sql) {
    if (interval != "none") stop("Database prediction does not support intervals.", call. = FALSE)
    return(predict_keyed_sql(object, newdata, conn, table, path, row_id, prediction_table))
  }
  if (interval != "none") validate_confidence_level(level)
  if (sample == "original") require_inference_rows(object)
  if (is.null(newdata) && isTRUE(object[["retain"]] %in% c("inference", "all"))) {
    if (interval != "none") stop("Retained-sample prediction requires interval = 'none'.", call. = FALSE)
    return(retained_row_values(object, "__dbreg_saved_fit", ".fitted", na.rm = sample == "estimation"))
  }
  strategy = object[["strategy"]]
  if (!is.null(newdata) && length(object[["fe"]]) && !inherits(object, "dbivreg") &&
      isTRUE(object[["retain"]] %in% c("prediction", "all"))) {
    if (interval != "none") stop("Learned-effect predictions require interval = 'none'.", call. = FALSE)
    if (is.null(object[["fitted_state"]])) {
      stop("The fit has no portable learned-effect state. Inspect prediction_diagnostics.", call. = FALSE)
    }
    return(predict_learned_effects(object, as.data.frame(newdata)))
  } 

  # Demean strategy doesn't support intervals (FE uncertainty not available)
  if (strategy == "demean" && interval != "none") {
    message(
      "Confidence/prediction intervals not supported for demean strategy ",
      "(fixed effect uncertainty unavailable). Returning point predictions."
    )
    interval = "none"
  }

  stored_cells = is.null(newdata)
  if (stored_cells && !is.null(object[["offset"]])) {
    stop("Offset models require newdata for prediction.", call. = FALSE)
  }
  if (is.null(newdata)) {
    if (strategy == "compress" && !is.null(object[["prediction_cells"]])) {
      newdata = object[["prediction_cells"]]
    } else if (strategy == "compress" && !is.null(object[["data"]])) {
      newdata = object[["data"]]
    } else {
      stop("newdata is required for predictions, as dbreg does not retain the original data.")
    }
  }

  coordinate_cells = strategy == "compress" && !is.null(object[["compressed_design"]][["bindings"]])
  if (strategy == "compress" && !stored_cells && !coordinate_cells) newdata = prepare_composite_prediction(object, newdata)
  prediction_offset = linear_prediction_offset(object, newdata)
  cohort = NULL
  if (strategy %in% c("demean", "mundlak")) {
    newdata = as.data.frame(newdata)
    cohort = prepare_prediction_cohort(object, newdata)
    newdata = cohort[["data"]]
    prediction_offset = prediction_offset[cohort[["valid"]]]
    if (!nrow(newdata)) {
      blank = rep(NA_real_, length(cohort[["valid"]]))
      return(if (interval == "none") blank else data.frame(fit = blank, lwr = blank, upr = blank))
    }
  }

  # Extract common components from object
  betas = coef(object, fe = TRUE)
  fml = object[["fml"]]
  fe = object[["fe"]]
  yvar = object[["yvar"]]
  intercept = "(Intercept)" %in% names(betas)

  # Ensure FE columns are factors
  if (is.null(object[["linear_design"]]) && strategy == "compress") for (f in fe) {
    newdata[[f]] = factor(newdata[[f]])
  }

  if (coordinate_cells) {
    mm = varying_compressed_prediction_matrix(object, newdata)
  } else if (!is.null(object[["linear_design"]])) {
    cells = object[["compressed_design"]]
    if (stored_cells && !is.null(cells)) {
      names(newdata) = cells[["data_names"]]
      newdata = linear_matrix_data(newdata)
      if (!is.null(cells[["columns"]])) {
        mm = as.matrix(newdata[, setdiff(cells[["columns"]], "(Intercept)"), drop = FALSE])
        if ("(Intercept)" %in% cells[["columns"]]) mm = cbind("(Intercept)" = 1, mm)
        mm = mm[, cells[["columns"]], drop = FALSE]
      } else {
        mm = sparse.model.matrix(cells[["formula"]], newdata, contrasts.arg = cells[["contrasts"]])[, cells[["positions"]], drop = FALSE]
      }
    } else {
      mm = linear_prediction_matrix(object, newdata)
    }
  } else {
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
  # (sparse.model.matrix() can't handle backticked, non-syntactic names)
  mm = tryCatch(
    sparse.model.matrix(fml_mm, data = newdata),
    error = function(e) model.matrix(fml_mm, data = newdata)
  )
  if (NROW(i_cols)) {
    is_i = colnames(mm) %in% paste0("`", i_cols[["name"]], "`")
    colnames(mm)[is_i] = gsub("`", "", colnames(mm)[is_i], fixed = TRUE)
  }
  if (strategy == "mundlak") {
    # Group means of each design column, one set per FE
    mm_means = mapply(
      function(term, f) weighted_group_mean(as.numeric(mm[, term]), newdata[[f]], prediction_weights(newdata, object[["weights"]])),
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

  }

  invalid_rows = attr(mm, "dbreg_invalid_rows", exact = TRUE)
  if (strategy == "demean") {
    has_y = yvar %in% names(newdata)
    working = as.matrix(mm)
    if (has_y) working = cbind(working, as.numeric(newdata[[yvar]]) - prediction_offset)
    projected = if (length(unlist(object[["varying_slopes"]]))) {
      demean_prediction_coordinates(working, prediction_fe_groups(object, newdata),
        lapply(object[["varying_slopes"]], function(vars) newdata[vars]), prediction_weights(newdata, object[["weights"]]))
    } else demean_prediction_data(working, prediction_fe_groups(object, newdata),
      prediction_weights(newdata, object[["weights"]]))
    mm = projected[, seq_len(ncol(mm)), drop = FALSE]
    if (has_y) {
      y_group_mean = as.numeric(newdata[[yvar]]) - prediction_offset - projected[, ncol(projected)]
    } else {
      y_group_mean = 0
      message(sprintf("Outcome '%s' not found in newdata. Returning within-group predictions (demean strategy only).", yvar))
    }
  }

  # Predict in the fitted moment coordinates, avoiding large
  # intercept/slope cancellation when the raw origins are large.
  prediction_coordinates = object[["prediction_coordinates"]]
  if (!is.null(prediction_coordinates)) {
    mm = center_compressed_design(mm, prediction_coordinates[["centers"]])
    betas = prediction_coordinates[["beta"]]
  }
  # Generate predictions
  fit = as.vector(mm %*% betas)
  if (!is.null(prediction_coordinates)) fit = fit + prediction_coordinates[["y_center"]]
  if (strategy == "demean") {
    fit = fit + y_group_mean
  }
 
  fit = fit + prediction_offset
  if (!is.null(invalid_rows)) fit[invalid_rows] = NA_real_
  if (interval != "none") {
    vcovm = vcov(object, fe = TRUE)
    if (!is.null(prediction_coordinates)) vcovm = prediction_coordinates[["covariance"]]
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
    t_crit = qt(1 - a, linear_test_df(object))
    lwr = fit - t_crit * ses
    upr = fit + t_crit * ses
    fit = data.frame(fit = fit, lwr = lwr, upr = upr)
  }

  if (!is.null(cohort)) {
    if (is.data.frame(fit)) {
      restored = as.data.frame(matrix(NA_real_, length(cohort[["valid"]]), ncol(fit)))
      names(restored) = names(fit); restored[cohort[["valid"]], ] = fit
    } else {
      restored = rep(NA_real_, length(cohort[["valid"]])); restored[cohort[["valid"]]] = fit
    }
    fit = restored
  }
  return(fit)
}

#' Variance-covariance matrix for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param fe Include covariance rows and columns for explicit fixed effects?
#'   Default is `FALSE`, matching [coef()].
#' @param vcov Post-fit covariance type or cluster formula for a retained OLS fit.
#' @param cluster One or two cluster terms from the original fitted source.
#' @param ssc Optional post-fit small-sample correction policy.
#' @param vcov_fix Repair an indefinite two-way covariance matrix.
#' @param ... Additional arguments (currently unused, except to capture
#'   superseded arguments).
#' @inherit coef.dbreg examples
#' @export
# Positions of the reported (non-FE) coefficients in a fitted-order covariance.
linear_reported_positions = function(object, n) {
  if (!length(object[["fe"]])) return(seq_len(n))
  if (!is.null(object[["linear_design"]])) return(visible_linear_positions(object, FALSE))
  match(names(coef(object)), rownames(object[["coeftable"]]))
}

vcov.dbreg = function(object, fe = FALSE, ..., vcov = NULL, cluster = NULL,
                      ssc = NULL, vcov_fix = TRUE) {
  dots = list(...)
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning('The `fes` argument has been superseded by `fe` (without the "s") and will be deprecated in a future `dbreg` release.\n')
    }
  }
  requested = !is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE)
  if (requested && is.null(vcov) && is.null(cluster)) {
    request = object[["postfit_request"]]
    vcov = request[["vcov"]]
    cluster = request[["cluster"]]
  }
  out = if (requested) retained_linear_covariance(object, vcov, cluster, ssc, vcov_fix) else object[["vcov"]]
  fields = c("type", "n_clusters", "n_params", "ssc", "df.K", "df.t", "rss", "tss", "repair")
  metadata = lapply(fields, function(name) attr(out, name))
  names(metadata) = fields
  # Public labels may differ from SQL aliases. Both objects share the fitted
  # parameter order; without FE filtering, retain positions even if labels repeat.
  dimnames(out) = rep(list(rownames(object[["coeftable"]])), 2L)
  visible = linear_reported_positions(object, nrow(out))
  # As in fixest, repair two-way covariance on the reported slope coefficients:
  # fitting coordinates can include absorbed FE directions that need not be PSD.
  repair = metadata[["repair"]]
  if (requested && !inherits(object, "dbglm") && identical(metadata[["type"]], "twoway") && is.null(repair)) {
    block = retained_vcov_psd(out[visible, visible, drop = FALSE], vcov_fix)
    out[visible, visible] = block
    repair = attr(block, "repair")
  }
  if (!isTRUE(fe) && length(object[["fe"]])) out = out[visible, visible, drop = FALSE]
  for (name in names(metadata)) if (!is.null(metadata[[name]])) attr(out, name) = metadata[[name]]
  if (!is.null(repair)) attr(out, "repair") = repair
  out
}

#' Confidence intervals for dbreg objects
#'
#' @param object A `dbreg` object.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level One finite numeric confidence level strictly between zero and
#'   one. Default is 0.95.
#' @param fe Should the fixed effects be included? Default is `FALSE`.
#' @param vcov,cluster,ssc,vcov_fix Optional post-fit covariance arguments
#'   forwarded to [summary.dbreg()].
#' @param ... Additional arguments. Currently unused, except to capture
#' superseded arguments.
#' @inherit coef.dbreg examples
#' @importFrom stats qt
#' @export
confint.dbreg = function(object, parm, level = 0.95, fe = FALSE, ...,
                         vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE) {
  validate_confidence_level(level)

  # superseded args handled through ...
  dots = list(...)
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning('The `fes` argument has been superseded by `fe` (without the "s") and will be deprecated in a future `dbreg` release.\n')

    }
  }

  object = summary.dbreg(object, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
  ct = visible_coefficient_table(object, fe)
  cf = ct[, "estimate"]
  ses = ct[, "std.error"]
  a = (1 - level) / 2
  critical = if (normal_inference(object)) stats::qnorm(1 - a) else qt(1 - a, linear_test_df(object))
  ci = cbind(cf - critical * ses, cf + critical * ses)
  rownames(ci) = rownames(ct)
  colnames(ci) = sprintf("%.1f %%", 100 * c(a, 1 - a))
  
  if (!missing(parm)) {
    ci = ci[parm, , drop = FALSE]
  }
  
  ci
}

validate_confidence_level = function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must lie strictly between zero and one.", call. = FALSE)
  }
  invisible(NULL)
}

