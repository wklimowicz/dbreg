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
      warning("The `fes` argument has been superseded by `fe` (without the \"s\") and will be deprecated in a future `dbreg` release.\n")
    }
  }
  if (isTRUE(fe) && !isTRUE(collin) && !is.null(object[["all_coefficients"]])) return(object[["all_coefficients"]])

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
#' @param conn,table,path Optional database connection and prediction source,
#'   following [dbreg()] source conventions. Supply exactly one of `newdata`,
#'   `table`, or `path`. DuckDB and native Databricks support database output.
#' @param row_id Explicit unique, nonmissing prediction key column names.
#'   Supplying database arguments or `row_id` selects keyed database output.
#' @param prediction_table Optional caller-named session table receiving output.
#'   Requires a caller-owned connection and must not already exist. Otherwise
#'   the result owns its temporary table; see [dbreg-prediction].
#' @param sample Without `newdata`, `"estimation"` returns the fitted rows and
#'   `"original"` includes omitted original rows as missing. Requires `retain="all"`;
#'   R inputs return ordered vectors and SQL inputs return keyed relations.
#' @section Database prediction:
#' Returns a `dbreg_prediction` relation containing the original typed keys and
#'   a prediction column. Its `prediction_column` field names that column.
#'   Rows have no implicit order. `head()` previews and `as.data.frame()` collects
#'   explicitly. Collection returns exact integer keys as text if the driver
#'   would otherwise round them; database columns retain their original types.
#'   Missing inputs and unknown/non-estimable FE combinations yield SQL `NULL`.
#'   Ordinary unseen factor levels raise an error. Learned effects require
#'   `retain = "prediction"`; database dictionaries require their owning session.
#'   Intervals, IV and Mundlak database output are currently unsupported.
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
  conn = NULL,
  table = NULL,
  path = NULL,
  row_id = NULL,
  prediction_table = NULL,
  sample = c("estimation", "original"),
  ...
) {
  check_unused_args(list(...))
  interval = match.arg(interval)
  sample = match.arg(sample)
  if (any(!vapply(list(conn, table, path, row_id, prediction_table), is.null, logical(1)))) {
    if (interval != "none") stop("Database prediction output does not support intervals.", call. = FALSE)
    return(predict_sql(object, newdata, conn, table, path, row_id, prediction_table))
  }
  if (is.null(newdata) && !is.null(object[["fitted_state"]][["inference"]])) {
    if (interval != "none") stop("Retained sample prediction currently supports interval='none'.", call. = FALSE)
    return(predict_retained_sample(object, sample = sample))
  }
  if (is.null(newdata) && sample == "original") stop("Original-sample prediction requires retain='all'.", call. = FALSE)
  if (is.null(newdata) && !is.null(object[["offset"]])) stop("Offset models require newdata or retain='all' for prediction.", call. = FALSE)
  stored_compression = identical(object[["strategy"]], "compress") &&
    (is.null(newdata) || identical(newdata, object[["data"]]))
  if (!is.null(object[["fitted_state"]][["prediction"]]) && !stored_compression) return(predict_linear_state(object, newdata, interval))
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

  required = unique(c(
    all.vars(stats::formula(fml, lhs = 0, rhs = 1)), object[["offset"]],
    unlist(object[["composite_fes"]]), setdiff(fe, names(object[["composite_fes"]])),
    unlist(object[["varying_slopes"]]),
    if (strategy %in% c("demean", "mundlak")) object[["weights"]],
    if (strategy == "demean" && yvar %in% names(newdata)) yvar
  ))
  if (stored) required = names(newdata)
  missing = setdiff(required, names(newdata))
  if (length(missing)) stop("newdata is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  n_new = nrow(newdata)
  valid = rep(TRUE, n_new)
  for (name in required) {
    value = newdata[[name]]
    valid = valid & !is.na(value)
    if (is.numeric(value) && !inherits(value, "integer64")) valid = valid & is.finite(value)
  }
  if (!any(valid)) {
    blank = rep(NA_real_, n_new)
    return(if (interval == "none") blank else data.frame(fit = blank, lwr = blank, upr = blank))
  }
  newdata = newdata[valid, , drop = FALSE]
  native_data = if (!is.null(object[["linear_design"]])) newdata else NULL
  newdata = prediction_factors(object, newdata)

  # Ensure FE columns are factors
  for (f in fe) {
    if (!is.factor(newdata[[f]])) newdata[[f]] = factor(newdata[[f]])
  }

  # Rebuild the design from the model's own terms (so factors, interactions and
  # no-intercept models expand exactly as in estimation), then align its
  # columns to the coefficients by name
  terms_rhs = attr(terms(formula(fml, lhs = 0, rhs = 1)), "term.labels")
  i_cols = object[["i_columns"]]
  if (strategy %in% c("compress", "moments")) {
    # FE enter as dummies (compress)
    prediction_trends = object[["compressed_design"]][["prediction_trends"]]
    if (length(prediction_trends)) newdata = compressed_trend_data(newdata, prediction_trends, object[["compressed_design"]][["trend_centers"]])
    nuisance = unlist(lapply(seq_along(fe), function(i) {
      slopes = object[["varying_slopes"]][[i]]
      if (length(prediction_trends)) slopes = names(prediction_trends)[match(slopes, prediction_trends)]
      if (length(slopes)) paste(vapply(fe[i], function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1)),
        vapply(slopes, function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1)),
        sep = ":"
      ) else character()
    }), use.names = FALSE)
    terms_rhs = c(terms_rhs, fe, nuisance)
  }
  if (!length(terms_rhs)) {
    terms_rhs = "1"
  }
  fml_mm = reformulate(terms_rhs, intercept = intercept)
  if (!is.null(object[["linear_design"]])) {
    mm = native_linear_matrix(object, native_data, stored = stored && identical(object[["strategy"]], "compress"))
  } else if (stored && !is.null(object[["compressed_design"]])) {
    spec = object[["compressed_design"]]
    data = expand_newdata_i(stats::formula(fml, lhs = 0, rhs = 1), newdata, i_cols, object[["i_specs"]])[["data"]]
    names(data) = unname(spec[["source_names"]][names(data)])
    data = compressed_trend_data(data, spec[["trend_sources"]], spec[["trend_centers"]])
    mm = fitted_model_matrix(spec[["formula"]], data)
    mm = mm[, names(spec[["columns"]]), drop = FALSE]
    colnames(mm) = unname(spec[["columns"]])
  } else {
    mm = build_newdata_design(fml_mm, newdata, object[["factor_levels"]], i_cols, object[["i_specs"]])
    for (alias in names(object[["compressed_design"]][["prediction_trends"]])) {
      colnames(mm) = gsub(alias, coef_var_name(object[["compressed_design"]][["prediction_trends"]][[alias]]), colnames(mm), fixed = TRUE)
    }
  }
  if (strategy == "mundlak") {
    # Group means of each design column, one set per FE
    mm_means = mapply(
      function(term, f) weighted_group_mean(
        as.numeric(mm[, term]), newdata[[f]],
        if (is.null(object[["weights"]])) NULL else newdata[[object[["weights"]]]]
      ),
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
  if (is.null(object[["linear_design"]])) {
    mm = if (length(betas)) mm[, names(betas), drop = FALSE] else matrix(0, nrow(newdata), 0L)
  }

  if (strategy == "demean") {
    # demean: sweep the FE out of the regressors (and recover the outcome's
    # group means, if available, to return level predictions)
    has_y = yvar %in% names(newdata)
    groups = lapply(fe, function(f) newdata[[f]])
    weights = if (is.null(object[["weights"]])) NULL else newdata[[object[["weights"]]]]
    slopes = lapply(object[["varying_slopes"]], function(vars) newdata[vars])
    sweep_fe = function(v) demean_prediction_data(v, groups, weights, slopes = slopes)
    mm = sweep_fe(mm)
    if (has_y) {
      y = as.numeric(newdata[[yvar]]) - linear_prediction_offset(object, newdata)
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
  fit = prediction_matrix_product(mm, betas)
  if (strategy == "demean") {
    fit = fit + y_group_mean
  }
  fit = fit + linear_prediction_offset(object, newdata)

  if (interval != "none") {
    vcovm = vcov(object, fe = TRUE)
    dof = attr(vcovm, "df.t")
    if (is.null(dof)) dof = object[["df_residual"]]
    if (interval == "confidence") {
      ses = sqrt(Matrix::rowSums((mm %*% vcovm) * mm))
    } else if (interval == "prediction") {
      # Use stored RSS from training to estimate sigma^2
      rss = attr(object[["vcov"]], "rss")
      if (is.null(rss)) {
        stop("Prediction intervals require RSS, which is not stored in this model.")
      }
      sig2 = rss / object[["df_residual"]]
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
#' @param ... Additional arguments are not supported.
#' @param vcov,cluster Optional post-fit covariance type or cluster formula/names.
#'   Requires `retain="all"`. Supported requests are IID, HC1, one-/two-way OLS
#'   and one-way GLM clustering. `^` combines typed columns within a cluster.
#' @param ssc Optional small-sample correction list; see [summary.dbreg()].
#' @param vcov_fix Repair nonpositive covariance eigenvalues when needed.
#' @details Without a request, returns the fitted covariance unchanged. Post-fit
#'   requests use saved scores on exactly the fitted sample. Missing cluster
#'   values error, and absent columns require refitting or a supported mapping.
#'   Ordinary post-fit requests only cover slope covariance (`fe=FALSE`).
#' @inherit coef.dbreg examples
#' @export
vcov.dbreg = function(object, fe = FALSE, ..., vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE) {
  check_unused_args(list(...))
  if (inherits(fe, "formula") || is.character(fe)) {
    if (!is.null(vcov)) stop("Specify one covariance request.", call. = FALSE)
    vcov = fe
    fe = FALSE
  }
  if (!is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE)) {
    if (!identical(fe, FALSE)) stop("Post-fit absorbed-effect covariance is not supported; use fe=FALSE.", call. = FALSE)
    return(inference_covariance(requested_inference_result(object,
      vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)))
  }
  out = object[["vcov"]]
  if (isTRUE(fe) && length(object[["fe"]]) && isTRUE(object[["postfit_covariance"]])) stop("Post-fit absorbed-effect covariance is not supported; use fe=FALSE.", call. = FALSE)
  if (!is.null(object[["inference_result"]])) {
    visible = if (!isTRUE(fe) && length(object[["fe"]])) names(coef(object, fe = fe)) else object[["inference_result"]][["coefficient_names"]]
    return(inference_covariance(object[["inference_result"]], visible))
  }
  fields = c("type", "n_clusters", "n_params", "rss", "tss", "df.K", "df.t", "ssc", "vcov_type", "cluster_var", "df_t", "repair")
  metadata = lapply(fields, function(nm) attr(out, nm))
  names(metadata) = fields
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
#' @param vcov,cluster,ssc,vcov_fix Optional post-fit covariance controls, as in
#'   [summary.dbreg()]. New requests require `retain="all"` and `fe=FALSE`.
#' @inherit coef.dbreg examples
#' @importFrom stats qt
#' @export
confint.dbreg = function(
  object, parm, level = 0.95, fe = FALSE,
  vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE, ...
) {
  # superseded args handled through ...
  dots = list(...)
  check_unused_args(dots, allowed = "fes")
  if (length(dots)) {
    if (!is.null(dots[["fes"]]) && !identical(fe, dots[["fes"]])) {
      fe = dots[["fes"]]
      warning("The `fes` argument has been superseded by `fe` (without the \"s\") and will be deprecated in a future `dbreg` release.\n")
    }
  }
  if (!is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE) || isTRUE(object[["postfit_covariance"]])) {
    if (!identical(fe, FALSE)) stop("Post-fit absorbed-effect intervals are not supported; use fe=FALSE.", call. = FALSE)
    V = vcov.dbreg(object, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
    return(inference_confint(object, V, if (missing(parm)) NULL else parm, level))
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
  if (!is.null(attr(object[["vcov"]], "df.t"))) df = attr(object[["vcov"]], "df.t")

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

#' Methods for database binomial, Poisson and Gaussian GLMs
#' @param object,x A `dbglm` model.
#' @param ... Additional method arguments are not supported and produce an error.
#'   New-row prediction with absorbed effects requires `retain = "prediction"`
#'   when fitting. Unknown or non-estimable FE combinations return `NA`.
#'   Database-backed dictionaries need their live owning session; see
#'   [close.dbreg()] for release and serialization behavior.
#' @name dbglm-methods
#' @importFrom stats logLik nobs df.residual predict
NULL

#' @rdname dbglm-methods
#' @param type For `logLik`, `"fixest"` reports the likelihood using the
#'   `fixest::feglm()` convention (the default); `"weighted"` returns the
#'   objective with unrounded observation weights. These differ for fractional
#'   binomial weights. Gaussian default likelihood is profile likelihood minus
#'   one, following R's family AIC convention; `"weighted"` omits that constant.
#'   For `predict`, `"link"` returns the linear predictor and
#'   `"response"` returns the response mean.
#' @export
logLik.dbglm = function(object, type = c("fixest", "weighted"), ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  value = if (type == "weighted" && !is.null(object[["weighted_loglik"]])) object[["weighted_loglik"]] else object[["loglik"]]
  structure(value,
    class = "logLik", df = if (isFALSE(object[["finite_mle"]])) NA_integer_ else object[["n_params_model"]],
    nobs = object[["nobs_orig"]], finite_mle = object[["finite_mle"]]
  )
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
#' @param vcov,cluster,ssc,vcov_fix Optional post-fit covariance controls, as in
#'   [summary.dbreg()]. New requests require `retain="all"`.
#' @export
confint.dbglm = function(
  object, parm, level = 0.95, fe = FALSE,
  vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE, ...
) {
  check_unused_args(list(...))
  if (!identical(fe, FALSE)) stop("Absorbed-effect confidence intervals are not supported; use `fe = FALSE`.", call. = FALSE)
  if (!is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE) || isTRUE(object[["postfit_covariance"]])) {
    V = vcov.dbreg(object, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
    return(inference_confint(object, V, if (missing(parm)) NULL else parm, level))
  }
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must lie strictly between zero and one.", call. = FALSE)
  }
  beta = stats::coef(object)
  if (missing(parm)) parm = names(beta)
  if (is.numeric(parm)) parm = names(beta)[parm]
  if (anyNA(parm) || any(!parm %in% names(beta))) stop("Unknown coefficient in `parm`.", call. = FALSE)
  probabilities = c((1 - level) / 2, (1 + level) / 2)
  critical = if (object[["family"]][["family"]] == "gaussian") stats::qt(probabilities, object[["df_t"]]) else stats::qnorm(probabilities)
  result = outer(sqrt(diag(stats::vcov(object)))[parm], critical) + beta[parm]
  colnames(result) = paste0(format(100 * c((1 - level) / 2, (1 + level) / 2)), " %")
  result
}

#' @rdname dbglm-methods
#' @param newdata Data frame containing regressors, fitted FE keys and varying
#'   slope variables, and any fitted offset column. Outcomes and fitting weights
#'   are not required. Missing inputs preserve their positions as `NA` predictions.
#' @param conn,table,path Optional database connection and prediction source, as
#'   in [predict.dbreg()]. Database output supports DuckDB and native Databricks.
#' @param row_id Unique, nonmissing prediction key columns. This or any database
#'   argument selects keyed database output; see [dbreg-prediction].
#' @param prediction_table Optional caller-named output session table. Requires
#'   a caller-owned connection and must not already exist.
#' @param sample Without `newdata`, returns fitted estimation rows or original
#'   rows with omitted positions missing; requires `retain="all"`.
#' @export
predict.dbglm = function(
  object, newdata = NULL, type = c("link", "response"),
  conn = NULL, table = NULL, path = NULL, row_id = NULL, prediction_table = NULL,
  sample = c("estimation", "original"), ...
) {
  check_unused_args(list(...))
  type = match.arg(type)
  sample = match.arg(sample)
  if (any(!vapply(list(conn, table, path, row_id, prediction_table), is.null, logical(1)))) {
    return(predict_sql(object, newdata, conn, table, path, row_id, prediction_table, type))
  }
  if (is.null(newdata) && !is.null(object[["fitted_state"]][["inference"]])) return(predict_retained_sample(object, type, sample))
  if (is.null(newdata)) stop("`newdata` is required; dbglm does not retain input rows.", call. = FALSE)
  if (length(object[["fe"]]) && is.null(object[["fitted_state"]][["prediction"]])) {
    stop("New-row prediction for absorbed GLMs requires fitting with retain='prediction'.", call. = FALSE)
  }
  newdata = as.data.frame(newdata)
  for (a in names(object[["factor_sources"]])) {
    source = object[["factor_sources"]][[a]]
    if (!source %in% names(newdata)) stop("newdata is missing columns: ", source, call. = FALSE)
    newdata[[a]] = newdata[[source]]
  }
  state = object[["fitted_state"]]
  if (is.null(state[["prediction"]])) state = NULL
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) fitted_resource_connection(object[["resource_state"]])
  required = unique(c(
    object[["model_spec"]][["source_vars"]], object[["offset"]],
    unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(object[["varying_slopes"]])
  ))
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
  lookup = NULL
  if (!is.null(state)) {
    lookup = if (database) prediction_database_lookup(state, object[["resource_state"]], data) else NULL
    groups = if (!database) prediction_match_groups(state, data) else NULL
    known = if (database) !is.na(lookup[["gamma"]]) else if (length(groups))
      Reduce(`&`, lapply(groups, function(group) !is.na(group))) else rep(TRUE, nrow(data))
    if (!any(known)) return(prediction)
    if (!all(known)) {
      data = data[known, , drop = FALSE]
      groups = lapply(groups, function(group) group[known])
      if (database) lookup = lookup[known, , drop = FALSE]
      positions = positions[known]
    }
  }
  design = as.matrix(evaluate_model_spec(object[["model_spec"]], data))
  colnames(design) = object[["model_spec"]][["columns"]][["name"]]
  if (object[["has_intercept"]]) design = cbind("(Intercept)" = 1, design)
  beta = stats::coef(object)
  names(beta) = object[["design_names"]]
  eta = if (!is.null(state)) predict_from_fitted_state(state, design, data, groups = groups, lookup = lookup) else
    prediction_matrix_product(design[, names(beta), drop = FALSE], beta)
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
  data.frame(
    nobs = x[["nobs_orig"]], nobs.input = x[["nobs_input"]], nobs.dropped = x[["nobs_dropped"]],
    logLik = x[["loglik"]], deviance = x[["deviance"]], AIC = stats::AIC(x), BIC = stats::BIC(x),
    df.residual = x[["df_residual"]], n.params = x[["n_params_model"]], family = x[["family"]][["family"]],
    link = x[["family"]][["link"]], vcov.type = x[["vcov_type"]], converged = x[["converged"]],
    finite.mle = x[["finite_mle"]], iterations = x[["iterations"]], score.max = x[["score_max"]]
  )
}

#' @rdname dbglm-methods
#' @export
print.dbglm = function(x, ...) {
  family = x[["family"]]
  cat("GLM estimation, family = ", family[["family"]], "(", family[["link"]],
    "), Dep. Var.: ", x[["yvar"]], "\n",
    sep = ""
  )
  cat("Observations:", prettyNum(x[["nobs_orig"]], big.mark = ","), "\n")
  if (length(x[["n_fe_levels"]])) cat(
    "Fixed-effects:",
    paste(paste0(names(x[["n_fe_levels"]]), ": ", prettyNum(x[["n_fe_levels"]], big.mark = ",")), collapse = "; "), "\n"
  )
  cat("Standard Errors:", format_vcov_label(x[["vcov_type"]], x[["n_clusters"]]), "\n")
  if (isFALSE(x[["finite_mle"]])) cat("Boundary probability approximation; no finite MLE established; inference unavailable.\n")
  table = as.matrix(x[["coeftable"]])
  normal = family[["family"]] %in% c("binomial", "poisson")
  colnames(table) = c("Estimate", "Std. Error", if (normal) "z value" else "t value", if (normal) "Pr(>|z|)" else "Pr(>|t|)")
  print_coeftable(table)
  cat(
    "Log-Likelihood:", format(x[["loglik"]], digits = 6),
    "  Deviance:", format(x[["deviance"]], digits = 6), "  BIC:", format(stats::BIC(x), digits = 6), "\n"
  )
  if (!is.null(x[["null_loglik"]]) && is.finite(x[["null_loglik"]]) && x[["null_loglik"]] != 0) {
    cat("Pseudo R2:", format(1 - x[["loglik"]] / x[["null_loglik"]], digits = 6), "\n")
  }
  cat("IRLS iterations:", x[["iterations"]], "\n")
  invisible(x)
}

# Portable prediction state contains values only. Observation-sized relations
# remain fit-owned and are released by the usual fitting resource manager.
prediction_state_limit = function() {
  limit = getOption("dbreg.prediction_max_entries", 2e6)
  validate_moment_limits(limit, limit)
  limit
}

prediction_state_overflow = function() {
  stop(structure(list(
    message = paste0(
      "Portable prediction state exceeds its transfer or matrix limit; ",
      "database-backed recovery is unavailable for this design or budget. Use retain='none' for a compact fit."
    ),
    call = NULL
  ), class = c("dbreg_prediction_overflow", "error", "condition")))
}

prediction_key_part = function(value, type) {
  if (type == "integer") {
    if (is.numeric(value) && !inherits(value, "integer64")) {
      bad = !is.na(value) & (!is.finite(value) | value != floor(value) | abs(value) > 2^53)
      if (any(bad)) stop("Integer FE keys must be exact; use integer64 or decimal strings for large keys.", call. = FALSE)
      value = format(value, scientific = FALSE, trim = TRUE, digits = 22)
    } else value = as.character(value)
    valid = !is.na(value) & grepl("^[+-]?[0-9]+$", value)
    value[!valid] = NA_character_
    value = sub("^\\+", "", value)
    value = sub("^(-?)0+([0-9])", "\\1\\2", value)
    value[value %in% c("-0", "+0")] = "0"
  } else if (type %in% c("numeric", "float")) {
    if (!is.numeric(value) || inherits(value, "integer64")) {
      stop("Numeric FE keys require numeric newdata columns.", call. = FALSE)
    }
    if (type == "float") value = readBin(writeBin(as.double(value), raw(), size = 4L), "double", n = length(value), size = 4L)
    value[value == 0 & !is.na(value)] = 0
    value = ifelse(is.finite(value), sprintf("%.17g", value), NA_character_)
  } else if (type == "logical") {
    if (!is.logical(value)) stop("Logical FE keys require logical newdata columns.", call. = FALSE)
    value = ifelse(is.na(value), NA_character_, ifelse(value, "TRUE", "FALSE"))
  } else if (type == "date") {
    if (!inherits(value, "Date") && !is.character(value)) stop("Date FE keys require Date or ISO date-string columns.", call. = FALSE)
    value = as.character(value)
  } else value = as.character(value)
  ifelse(is.na(value), NA_character_, paste0(nchar(enc2utf8(value), type = "bytes"), ":", enc2utf8(value)))
}

prediction_tuple_key = function(data, types) {
  parts = Map(prediction_key_part, data, types)
  key = do.call(paste0, unname(parts))
  missing = Reduce(`|`, lapply(parts, is.na))
  key[missing] = NA_character_
  key
}

prediction_key_fields = function(conn, components) {
  components = lapply(components, unique)
  columns = lapply(seq_along(components), function(i) paste0("__dbglm_key_fe", i, "_", seq_along(components[[i]])))
  fields = unlist(lapply(seq_along(components), function(i) {
    paste(as.character(DBI::dbQuoteIdentifier(conn, components[[i]])), "AS", columns[[i]])
  }), use.names = FALSE)
  list(fields = fields, columns = columns)
}

prediction_key_type = function(type) {
  type = tolower(type)
  if (grepl("^(tinyint|smallint|int|integer|bigint|hugeint|u(tinyint|smallint|integer|bigint|hugeint)|decimal\\([0-9]+,0\\))$", type)) "integer" else if (type %in% c("double", "double precision")) "numeric" else if (type %in% c("float", "real")) "float" else if (type %in% c("boolean", "bool")) "logical" else if (grepl("^(varchar|char|string|enum)", type)) "character" else if (type == "date") "date" else NA_character_
}

prediction_key_metadata = function(conn, retained, components) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  sql_types = backend_get_query(conn, paste0(
    "SELECT ", paste0("typeof(", q(components), ") AS t", seq_along(components), collapse = ","),
    " FROM ", retained, " LIMIT 1"
  ))
  types = vapply(sql_types, function(type) prediction_key_type(type[[1L]]), character(1))
  if (anyNA(types)) {
    unknown = tolower(sql_types[[which(is.na(types))[[1L]]]][[1L]])
    stop("Portable FE prediction does not support key SQL type '", unknown, "'.", call. = FALSE)
  }
  list(types = unname(types), sql_types = vapply(sql_types, function(type) type[[1]], character(1)))
}

prediction_fe_dictionary = function(conn, retained, fe, components, limit, public_components = components) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  metadata = prediction_key_metadata(conn, retained, components)
  types = metadata[["types"]]
  fields = vapply(seq_along(components), function(i) {
    value = q(components[i])
    if (types[i] == "integer") value = paste0("CAST(", value, " AS STRING)")
    paste(value, "AS", paste0("k", i))
  }, character(1))
  max_rows = floor(limit / (length(components) + 1L))
  if (max_rows < 1L) prediction_state_overflow()
  dictionary = backend_get_query(conn, paste0(
    "SELECT DISTINCT ", fe, " AS id,", paste(fields, collapse = ","),
    " FROM ", retained, " ORDER BY id LIMIT ", sprintf("%.0f", max_rows + 1)
  ))
  if (nrow(dictionary) > max_rows) prediction_state_overflow()
  keys = prediction_tuple_key(dictionary[paste0("k", seq_along(components))], types)
  if (anyNA(keys) || anyDuplicated(keys) || anyDuplicated(dictionary[["id"]])) {
    stop("Could not preserve the fitted FE key dictionary exactly.", call. = FALSE)
  }
  list(
    columns = public_components, types = unname(types), key = keys, id = dictionary[["id"]],
    sql_types = unname(metadata[["sql_types"]])
  )
}

dbglm_prediction_state = function(...) build_prediction_state(...)

build_prediction_state = function(conn, working, retained, fe, components, x, design_names,
                                  center, trend_aliases, sizes, create, drop, n,
                                  dictionary_components = components, beta = NULL,
                                  expected_rss = NULL, constant = 0, fitted_beta = beta,
                                  recovery_control = list(tol = 1e-10, max_iter = 100L)) {
  tryCatch(
    build_portable_prediction_state(
      conn, working, retained, fe, components, x, design_names,
      center, trend_aliases, sizes, create, drop, n, dictionary_components, beta, expected_rss, constant
    ),
    dbreg_prediction_overflow = function(error) {
      if (!length(fe) || length(fe) > 2L || length(unlist(attr(fe, "slopes"))) || is.null(fitted_beta)) stop(error)
      if (length(fe) == 2L) {
        return(build_database_two_fe_state(
          conn, working, retained, fe, components, x, design_names,
          center, create, drop, n, dictionary_components, fitted_beta, expected_rss, constant, recovery_control
        ))
      }
      build_database_prediction_state(
        conn, working, retained, fe, components, x, design_names,
        center, create, n, dictionary_components, fitted_beta, expected_rss, constant
      )
    }
  )
}

build_portable_prediction_state = function(conn, working, retained, fe, components, x, design_names,
                                           center, trend_aliases, sizes, create, drop, n,
                                           dictionary_components = components, beta = NULL,
                                           expected_rss = NULL, constant = 0, precomputed = NULL) {
  backend_info = detect_backend(conn)
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  limit = prediction_state_limit()
  limits = projection_limits()
  max_rows = limits[["summary_rows"]]
  validate_moment_limits(max_rows, limit)
  used = 0
  bounded_query = function(sql, width) {
    remaining = min(max_rows - used, floor(limit / width))
    if (remaining < 1) prediction_state_overflow()
    result = backend_get_query(conn, paste0(sql, " LIMIT ", sprintf("%.0f", remaining + 1)))
    if (nrow(result) > remaining) prediction_state_overflow()
    used = used + nrow(result)
    result
  }
  ops = lapply(seq_along(fe), function(i) list(f = fe[i], v = "1.0", group = i, variable = NULL, center = 0))
  for (i in seq_along(fe)) for (variable in attr(fe, "slopes")[[i]]) {
    name = names(trend_aliases)[match(variable, trend_aliases)]
    centers = attr(trend_aliases, "centers")
    ops[[length(ops) + 1L]] = list(
      f = fe[i], v = variable, group = i,
      variable = name, center = if (is.null(centers)) 0 else centers[[name]]
    )
  }
  source = if (is.null(precomputed)) paste0("(", backend_sql(conn, paste0(
    "SELECT * EXCLUDE(w), __dbglm_a AS w,",
    "(eta-__dbglm_o) AS __dbglm_prediction_target FROM ", working
  )), ")") else NULL
  columns = c(x, "__dbglm_prediction_target")
  op_sizes = sizes[vapply(ops, `[[`, integer(1), "group")]
  total = sum(op_sizes) + length(x)
  rest = total - max(op_sizes)
  if (!fe_dense_feasible(c(op_sizes, rep(1, length(x))), length(columns), limits) || total * max(1, rest, length(columns)) > limit) {
    prediction_state_overflow()
  }
  if (!is.null(precomputed)) {
    info = precomputed[["info"]]
    sys = precomputed[["sys"]]
    B = precomputed[["B"]]
  } else if (length(fe) >= 2L && length(ops) == length(fe)) {
    payload = sum(sizes) * length(columns)
    # Recovery's transient sufficient statistics use the shared transfer budget;
    # prediction_max_entries bounds the retained matrices and key dictionaries.
    row_cap = min(max_rows, floor((limits[["transfer_cells"]] - payload) / 5L) - 1L)
    if (row_cap < 1L) prediction_state_overflow()
    sql = fe_moments_sql(conn, source, fe, fe_value_terms(conn, columns, backend_info = backend_info),
      limit = row_cap + 1, compact = TRUE
    )
    res = read_fe_moments(conn, sql)
    if (nrow(res) > row_cap) prediction_state_overflow()
    used = nrow(res)
    moments = unpack_fe_moments(res, length(fe), length(columns), 0L)
    if (is.null(moments)) prediction_state_overflow()
    info = moments[["info"]]
    sys = moments[["sys"]]
    B = moments[["B"]]
  } else {
    info = lapply(ops, function(op) bounded_query(paste0(
      "SELECT ", op[["f"]], " AS id,",
      sum_sql(paste0("w*(", op[["v"]], ")*(", op[["v"]], ")")), " AS h,",
      paste(paste(sum_sql(sprintf("w*(%s)*(%s)", op[["v"]], columns)), "AS", paste0("b", seq_along(columns))), collapse = ","),
      " FROM ", source, " GROUP BY ", op[["f"]], " ORDER BY ", op[["f"]]
    ), length(columns) + 2L))
    pairs = list()
    for (i in seq_along(ops)) if (i > 1L) for (j in seq_len(i - 1L)) {
      oi = ops[[i]]
      oj = ops[[j]]
      edge = bounded_query(paste0(
        "SELECT ", oi[["f"]], " AS i,", oj[["f"]], " AS j,",
        sum_sql(paste0("w*(", oi[["v"]], ")*(", oj[["v"]], ")")), " AS h FROM ", source,
        " GROUP BY ", paste(unique(c(oi[["f"]], oj[["f"]])), collapse = ",")
      ), 3L)
      pairs[[length(pairs) + 1L]] = list(
        i = i, j = j, h = as.numeric(edge[["h"]]),
        ir = match(edge[["i"]], info[[i]][["id"]]), jr = match(edge[["j"]], info[[j]][["id"]])
      )
    }
    sys = fe_system(info, pairs)
    if (is.null(sys)) prediction_state_overflow()
    B = do.call(rbind, lapply(info, function(group) as.matrix(group[paste0("b", seq_along(columns))])))
  }
  cross = if (!is.null(precomputed)) precomputed[["cross"]] else {
    global = backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns, "w"), collapse = ","), "FROM", source))
    dbglm_read_information(global, length(columns))
  }
  p = length(x)
  fe_n = sys[["n"]]
  fe_system = sys
  if (p) {
    old_rest = sys[["rest"]]
    bx = B[, seq_len(p), drop = FALSE]
    sys[["Cm"]] = cbind(sys[["Cm"]], bx[sys[["big"]], , drop = FALSE])
    sys[["E"]] = rbind(
      cbind(sys[["E"]], bx[old_rest, , drop = FALSE]),
      cbind(t(bx[old_rest, , drop = FALSE]), cross[seq_len(p), seq_len(p), drop = FALSE])
    )
    sys[["rest"]] = c(old_rest, fe_n + seq_len(p))
    sys[["n"]] = fe_n + p
    sys[["diag"]] = c(sys[["diag"]], diag(cross)[seq_len(p)])
  }
  rhs = matrix(c(B[, p + 1L], if (p) cross[seq_len(p), p + 1L]), ncol = 1L)
  solver = fe_solver(sys, new.env(parent = emptyenv()), retain_null = TRUE)
  if (is.null(beta)) {
    rhs = solver[["project"]](rhs)
    coefficients = solver[["solve"]](rhs)
    # Recover the fitted GLM predictor, independently of its fitting projector.
    for (i in seq_len(3L)) coefficients = coefficients + solver[["solve"]](rhs - fe_multiply(sys, coefficients))
  } else {
    # OLS already supplies its fitted slopes. Solve only the nuisance equations
    # for y-Xb, preserving those slopes and their original estimation strategy.
    if (length(beta) != p || any(!is.finite(beta)) || is.null(expected_rss)) {
      stop("OLS prediction recovery requires fitted slopes and residual sum of squares.", call. = FALSE)
    }
    if (!is.null(precomputed[["alpha"]])) {
      alpha = precomputed[["alpha"]]
      gamma = as.numeric(alpha[, p + 1L] - alpha[, seq_len(p), drop = FALSE] %*% beta)
    } else {
      fe_solver = fe_solver(fe_system, new.env(parent = emptyenv()))
      residual_rhs = matrix(B[, p + 1L] - B[, seq_len(p), drop = FALSE] %*% beta, ncol = 1L)
      residual_rhs = fe_solver[["project"]](residual_rhs)
      gamma = fe_solver[["solve"]](residual_rhs)
      for (i in seq_len(3L)) gamma = gamma + fe_solver[["solve"]](residual_rhs - fe_multiply(fe_system, gamma))
    }
    coefficients = matrix(c(gamma, beta), ncol = 1L)
  }
  if (any(!is.finite(coefficients))) stop("Prediction coefficient recovery is nonfinite.", call. = FALSE)
  if (!is.null(precomputed)) {
    contrast = c(-beta, 1)
    recovered = as.numeric(crossprod(contrast, precomputed[["within"]] %*% contrast))
    recovery_error = abs(recovered - expected_rss)
    if (!is.finite(recovery_error) || recovery_error > 1e-6 * max(1, abs(expected_rss), cross[p + 1L, p + 1L])) {
      stop("Recovered OLS effects do not reproduce the fitted residual sum of squares.", call. = FALSE)
    }
  } else {
    frame = dbglm_coefficient_frame(info, coefficients[seq_len(fe_n), , drop = FALSE])
    lookup = dbglm_coefficient_lookup(conn, frame, create)
    ctes = dbglm_coefficient_ctes(lookup[["source"]], length(ops), "pred")
    joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN pred_g", i, " ON s.", ops[[i]][["f"]], "=pred_g", i, ".__glm_group"), character(1))
    terms = c(
      if (p) sprintf("(%s)*s.%s", dbglm_number(coefficients[fe_n + seq_len(p)]), x),
      vapply(seq_along(ops), function(i) sprintf("(%s)*pred_g%d.c1", if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]]), i), character(1))
    )
    residual = paste0("((", paste(terms, collapse = "+"), ")-s.__dbglm_prediction_target)")
    diagnostics = if (is.null(expected_rss)) paste0(
      "MAX(ABS(", residual,
      ")) AS error,MAX(ABS(s.__dbglm_prediction_target)) AS scale"
    ) else
      paste(sum_sql(paste0("s.w*POWER(", residual, ",2)")), "AS rss")
    error = backend_get_query(conn, paste0(
      "WITH ", paste(ctes, collapse = ","), " SELECT ", diagnostics, " FROM ",
      source, " s", paste(joins, collapse = "")
    ))
    for (table in lookup[["tables"]]) drop(table)
    if (is.null(expected_rss)) {
      recovery_error = error[["error"]]
      if (!is.finite(recovery_error) || recovery_error > 1e-7 * max(1, error[["scale"]])) {
        stop("Recovered prediction state does not reproduce the fitted linear predictor.", call. = FALSE)
      }
    } else {
      recovery_error = abs(error[["rss"]] - expected_rss)
      if (!is.finite(recovery_error) || recovery_error > 1e-6 * max(1, abs(expected_rss), cross[p + 1L, p + 1L])) {
        stop("Recovered OLS effects do not reproduce the fitted residual sum of squares (recovered=",
          signif(error[["rss"]], 9), ", fitted=", signif(expected_rss, 9), ").",
          call. = FALSE
        )
      }
    }
  }
  dictionaries = if (!is.null(precomputed)) precomputed[["dictionaries"]] else lapply(seq_along(fe), function(i) {
    prediction_fe_dictionary(conn, retained, fe[i], unique(dictionary_components[[i]]),
      floor(limit / length(fe)),
      public_components = unique(components[[i]])
    )
  })
  list(
    version = 1L, specification = list(
      design_names = design_names, design_center = unname(center),
      dictionaries = dictionaries, operations = lapply(ops, function(op) op[c("group", "variable", "center")])
    ),
    coefficients = list(
      joint = as.numeric(coefficients), starts = sys[["starts"]],
      ids = lapply(info, `[[`, "id"), fe_n = fe_n
    ),
    sample = list(nobs = n), inference = NULL, resources = NULL,
    prediction = list(
      null_basis = solver[["null_basis"]], tolerance = 1e-7,
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant
    )
  )
}

prediction_match_groups = function(state, data) {
  if (!identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
    match(prediction_tuple_key(data[dictionary[["columns"]]], dictionary[["types"]]), dictionary[["key"]])
  })
}

predict_from_fitted_state = function(state, design, data, groups = NULL, lookup = NULL) {
  if (!identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  if (identical(state[["prediction"]][["storage"]], "database")) {
    return(predict_from_database_state(state, design, lookup))
  }
  specification = state[["specification"]]
  coef = state[["coefficients"]]
  x = sweep(design[, specification[["design_names"]], drop = FALSE], 2L, specification[["design_center"]], "-")
  if (is.null(groups)) groups = prediction_match_groups(state, data)
  ops = specification[["operations"]]
  indices = values = vector("list", length(ops))
  for (i in seq_along(ops)) {
    group = ops[[i]][["group"]]
    dictionary = specification[["dictionaries"]][[group]]
    indices[[i]] = coef[["starts"]][i] + match(dictionary[["id"]][groups[[group]]], coef[["ids"]][[i]])
    values[[i]] = if (is.null(ops[[i]][["variable"]])) rep(1, nrow(data)) else
      as.numeric(data[[ops[[i]][["variable"]]]]) - if (is.null(ops[[i]][["center"]])) 0 else ops[[i]][["center"]]
  }
  prediction = prediction_matrix_product(x, coef[["joint"]][coef[["fe_n"]] + seq_len(ncol(x))])
  if (!is.null(state[["prediction"]][["constant"]])) prediction = prediction + state[["prediction"]][["constant"]]
  norm = rowSums(x * x)
  for (i in seq_along(ops)) {
    prediction = prediction + values[[i]] * coef[["joint"]][indices[[i]]]
    norm = norm + values[[i]]^2
  }
  null = state[["prediction"]][["null_basis"]]
  if (!is.null(null)) {
    chunk = max(1, min(4096, floor(prediction_state_limit() / ncol(null))))
    for (start in seq.int(1L, nrow(data), by = chunk)) {
      rows = seq.int(start, min(nrow(data), start + chunk - 1L))
      score = x[rows, , drop = FALSE] %*% null[coef[["fe_n"]] + seq_len(ncol(x)), , drop = FALSE]
      for (i in seq_along(ops)) score = score + values[[i]][rows] * null[indices[[i]][rows], , drop = FALSE]
      non_estimable = !is.finite(rowSums(score)) | sqrt(rowSums(score^2)) > state[["prediction"]][["tolerance"]] * sqrt(norm[rows])
      prediction[rows[non_estimable]] = NA_real_
    }
  }
  prediction[!is.finite(prediction)] = NA_real_
  prediction
}

# Use the original linear fit's slopes and RSS; only nuisance coefficients are
# recovered. Both model classes share the bounded FE kernel and prediction state.
finalize_linear_prediction = function(result, inputs) {
  if (!retain_prediction(inputs[["retain"]]) || !length(inputs[["fe"]])) return(result)
  reused = tryCatch(finalize_fwl_prediction(result, inputs),
    dbreg_prediction_overflow = function(error) NULL
  )
  if (!is.null(reused)) return(reused)
  conn = inputs[["conn"]]
  components = lapply(inputs[["prediction_fe_components"]], unique)
  from = inputs[["from_statement"]]
  source = paste0("(SELECT * ", from, ")")
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  design = sql_model_matrix_internal(inputs[["fml"]], conn, sub("^FROM\\s+", "", from, ignore.case = TRUE),
    expand = "all", fe_vars = inputs[["fe"]], intercept = inputs[["intercept"]], factor_refs = inputs[["factor_refs"]],
    reserved_names = sql_design_reserved(inputs)
  )
  register_factor_display(inputs, design[["factor_levels"]], design[["level_aliases"]], design[["col_names"]], design[["aliases"]], factor_classes = design[["factor_classes"]])
  map = coefficient_name_map(inputs[["name_map"]], inputs[["source_map"]])
  design_names = standardize_coef_names(restore_names(design[["aliases"]], map))
  beta = coef.dbreg(result)
  if (length(setdiff(names(beta), design_names))) stop("Could not align fitted OLS prediction slopes.", call. = FALSE)
  full_beta = stats::setNames(rep(0, length(design_names)), design_names)
  full_beta[names(beta)] = beta
  taken = unique(c(
    inputs[["yvar"]], inputs[["xvars"]], inputs[["fe"]], unlist(components),
    inputs[["weights"]], inputs[["cluster_var"]], inputs[["row_id_safe"]]
  ))
  fe = safe_identifiers(paste0("__dbglm_fe", seq_along(components)), taken = taken, fold_case = folds_identifier_case(conn))
  trend_sources = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  trend_columns = if (length(trend_sources)) safe_identifiers(
    paste0("__dbglm_prediction_t", seq_along(trend_sources)),
    taken = c(taken, fe),
    fold_case = folds_identifier_case(conn)
  ) else character()
  trend_aliases = stats::setNames(trend_columns, restore_names(trend_sources, inputs[["source_map"]]))
  attr(fe, "slopes") = lapply(inputs[["fe_slopes"]], function(vars) trend_columns[match(vars, trend_sources)])
  ranked = paste0("(", sql_rank_fixed_effects(conn, source, components, fe), ")")
  x = if (length(design_names)) paste0("__dbglm_prediction_x", seq_along(design_names)) else character()
  weights = if (is.null(inputs[["weights"]])) "1.0" else q(inputs[["weights"]])
  summary_terms = c(
    "COUNT(*) AS n", paste0("MAX(", weights, ") AS weight_max"),
    paste0("AVG(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)) AS y_center"),
    sprintf("AVG(CAST((%s) AS DOUBLE)) AS center%d", design[["select_exprs"]], seq_along(x)),
    sprintf("COUNT(DISTINCT %s) AS size%d", fe, seq_along(fe)),
    if (length(trend_sources)) sprintf(
      "MIN(%1$s)/2.0+MAX(%1$s)/2.0 AS trend_center%2$d",
      sql_cast_double(q(trend_sources), conn), seq_along(trend_sources)
    )
  )
  row = backend_get_query(conn, paste("SELECT", paste(summary_terms, collapse = ","), "FROM", ranked))
  n = as.numeric(row[["n"]])
  if (n != result[["nobs_orig"]]) stop("The OLS sample changed before prediction recovery.", call. = FALSE)
  center = vapply(seq_along(x), function(i) as.numeric(row[[paste0("center", i)]]), numeric(1))
  y_center = as.numeric(row[["y_center"]])
  if (any(!is.finite(c(center, y_center)))) stop("OLS prediction centering is nonfinite.", call. = FALSE)
  sizes = vapply(seq_along(fe), function(i) as.numeric(row[[paste0("size", i)]]), numeric(1))
  weight_max = as.numeric(row[["weight_max"]])
  weight_scale = if (weight_max < 1) 2^floor(log2(weight_max)) else 1
  trend_centers = vapply(seq_along(trend_sources), function(i) as.numeric(row[[paste0("trend_center", i)]]), numeric(1))
  attr(trend_aliases, "centers") = stats::setNames(trend_centers, names(trend_aliases))
  keys = prediction_key_fields(conn, components)
  fields = c(
    fe, keys[["fields"]],
    paste0("(CAST(", q(inputs[["yvar"]]), " AS DOUBLE)-(", dbglm_number(y_center), ")) AS eta"),
    "CAST(0.0 AS DOUBLE) AS __dbglm_o",
    paste0("CAST(", weights, " AS DOUBLE)/(", dbglm_number(weight_scale), ") AS __dbglm_a"),
    "CAST(1.0 AS DOUBLE) AS w",
    if (length(trend_sources)) paste0(
      "(", sql_cast_double(q(trend_sources), conn), "-(",
      dbglm_number(trend_centers), ")) AS ", trend_columns
    ),
    sprintf("(CAST((%s) AS DOUBLE)-(%s)) AS %s", design[["select_exprs"]], dbglm_number(center), x)
  )
  frame = paste0("(", backend_sql(conn, paste("SELECT", paste(fields, collapse = ","), "FROM", ranked)), ")")
  manager = fit_table_manager(conn, inputs[["sample_resources"]])
  public_components = lapply(components, restore_names, map = inputs[["source_map"]])
  result[["fitted_state"]] = build_prediction_state(conn, frame, frame, fe, public_components,
    x, design_names, center, trend_aliases, sizes, manager[["create"]], manager[["drop"]], n,
    dictionary_components = keys[["columns"]], beta = unname(full_beta),
    expected_rss = attr(result[["vcov"]], "rss") / weight_scale, constant = y_center,
    recovery_control = list(tol = inputs[["fe_tol"]], max_iter = inputs[["fe_max_iter"]])
  )
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_term_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  result
}

# Native typed IDs already arrived in compact FE summaries. Reuse them as
# prediction dictionaries, without scanning rows to rank the same keys again.
fwl_prediction_dictionary = function(ids, sql_type, column) {
  type = prediction_key_type(sql_type)
  if (is.na(type)) return(NULL)
  values = switch(type,
    numeric = as.numeric(ids),
    float = as.numeric(ids),
    logical = ifelse(toupper(ids) == "TRUE", TRUE, ifelse(toupper(ids) == "FALSE", FALSE, NA)),
    ids
  )
  keys = prediction_tuple_key(stats::setNames(list(values), column), type)
  if (anyNA(keys) || anyDuplicated(keys)) return(NULL)
  list(columns = column, types = type, key = keys, id = seq_along(ids), sql_types = sql_type)
}

finalize_fwl_prediction = function(result, inputs) {
  moments = inputs[["linear_prediction_moments"]]
  if (is.null(moments) || moments[["nobs"]] != result[["nobs_orig"]]) return(NULL)
  x = moments[["x"]]
  p = length(x)
  map = coefficient_name_map(inputs[["name_map"]], inputs[["source_map"]])
  design_names = standardize_coef_names(restore_names(x, map))
  fitted = coef.dbreg(result)
  if (length(setdiff(names(fitted), design_names))) return(NULL)
  beta = stats::setNames(rep(0, p), design_names)
  beta[names(fitted)] = fitted
  alpha = moments[["alpha"]]
  gamma = as.numeric(alpha[, 1L] - alpha[, seq_len(p) + 1L, drop = FALSE] %*% unname(beta))
  starts = cumsum(c(0, moments[["sizes"]]))
  nuisance_bound = sum(vapply(seq_along(moments[["sizes"]]), function(i)
    max(abs(gamma[starts[[i]] + seq_len(moments[["sizes"]][[i]])])), numeric(1)))
  x_bound = sqrt(pmax(diag(moments[["cross"]])[-1L], 0) / moments[["min_weight"]])
  prediction_bound = nuisance_bound + sum(abs(beta) * x_bound)
  response_scale = max(1, sqrt(max(moments[["cross"]][1L, 1L], 0) / moments[["sum_w"]]))
  # Keep the existing row certificate for numerically risky cancellation.
  if (!is.finite(prediction_bound) || prediction_bound > 1e8 * response_scale) return(NULL)
  components = lapply(inputs[["prediction_fe_components"]], unique)
  public = lapply(components, restore_names, map = inputs[["source_map"]])
  conn = inputs[["conn"]]
  source = paste0("(SELECT * ", inputs[["from_statement"]], ")")
  dictionaries = vector("list", length(components))
  info = moments[["info"]]
  for (i in seq_along(components)) {
    ids = info[[i]][["id"]]
    dictionary = if (length(components[[i]]) == 1L) {
      fwl_prediction_dictionary(ids, moments[["types"]][[i]], public[[i]])
    } else {
      key = paste0("CAST(", DBI::dbQuoteIdentifier(conn, inputs[["fe"]][[i]]), " AS VARCHAR)")
      dictionary = prediction_fe_dictionary(conn, source, key, components[[i]],
        floor(prediction_state_limit() / length(components)),
        public_components = public[[i]]
      )
      dictionary[["id"]] = match(as.character(dictionary[["id"]]), ids)
      if (length(dictionary[["id"]]) != length(ids) || anyNA(dictionary[["id"]])) return(NULL)
      dictionary
    }
    if (is.null(dictionary)) return(NULL)
    dictionaries[[i]] = dictionary
    info[[i]][["id"]] = seq_along(ids)
  }
  order = c(seq_len(p) + 1L, 1L)
  precomputed = list(
    info = info, sys = moments[["sys"]],
    B = moments[["B"]][, order, drop = FALSE], alpha = alpha[, order, drop = FALSE],
    cross = moments[["cross"]][order, order, drop = FALSE],
    within = moments[["within"]][order, order, drop = FALSE], dictionaries = dictionaries
  )
  state = tryCatch(
    build_portable_prediction_state(conn, NULL, NULL, inputs[["fe"]], public,
      x, design_names, moments[["centers"]][-1L], character(), moments[["sizes"]], NULL, NULL,
      moments[["nobs"]],
      beta = unname(beta),
      expected_rss = attr(result[["vcov"]], "rss") / moments[["weight_scale"]],
      constant = moments[["centers"]][[1L]], precomputed = precomputed
    ),
    dbreg_prediction_overflow = function(error) NULL
  )
  if (is.null(state)) return(NULL)
  result[["fitted_state"]] = state
  result
}

predict_linear_state = function(object, newdata, interval) {
  if (is.null(newdata)) stop("newdata is required; prediction retention does not keep input rows.", call. = FALSE)
  if (interval != "none") stop("Intervals with learned effects are not implemented; use interval='none'.", call. = FALSE)
  data = as.data.frame(newdata)
  rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
  state = object[["fitted_state"]]
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) fitted_resource_connection(object[["resource_state"]])
  required = unique(c(
    all.vars(rhs), object[["offset"]], unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable"))
  ))
  absent = setdiff(required, names(data))
  if (length(absent)) stop("newdata is missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  valid = rep(TRUE, nrow(data))
  for (name in required) {
    valid = valid & !is.na(data[[name]])
    if (is.numeric(data[[name]]) && !inherits(data[[name]], "integer64")) valid = valid & is.finite(data[[name]])
  }
  prediction = rep(NA_real_, nrow(data))
  if (!any(valid)) return(prediction)
  rows = data[valid, , drop = FALSE]
  positions = which(valid)
  lookup = if (database) prediction_database_lookup(state, object[["resource_state"]], rows) else NULL
  groups = if (!database) prediction_match_groups(state, rows) else NULL
  known = if (database) !is.na(lookup[["gamma"]]) else if (length(groups))
    Reduce(`&`, lapply(groups, function(group) !is.na(group))) else rep(TRUE, nrow(rows))
  if (!any(known)) return(prediction)
  if (!all(known)) {
    rows = rows[known, , drop = FALSE]
    groups = lapply(groups, function(group) group[known])
    if (database) lookup = lookup[known, , drop = FALSE]
    positions = positions[known]
  }
  levels = object[["factor_levels"]]
  levels = levels[intersect(names(levels), c(all.vars(rhs), names(explicit_factor_sources(attr(stats::terms(rhs), "term.labels")))))]
  design = as.matrix(build_newdata_design(rhs, rows, levels, object[["i_columns"]], object[["i_specs"]]))
  absent = setdiff(state[["specification"]][["design_names"]], colnames(design))
  if (length(absent)) stop("Could not reconstruct fitted OLS design: ", paste(absent, collapse = ", "), call. = FALSE)
  prediction[positions] = predict_from_fitted_state(state, design, rows, groups = groups, lookup = lookup) + linear_prediction_offset(object, rows)
  prediction
}

# SQL projections recover learned effects and regressor relations without
# collecting FE levels. Reusing native table names bounds the number of live
# intermediates; the fit manager owns all allocations until final adoption.
prediction_two_fe_projection = function(conn, source, fe, targets, create, drop, tol = 1e-10, max_iter = 100L) {
  backend_info = detect_backend(conn)
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  native = isTRUE(backend_info[["is_databricks"]])
  names = paste0("r", seq_along(targets))
  scale_row = backend_get_query(conn, paste("SELECT", paste(sprintf("MAX(ABS(%s)) AS s%d", targets, seq_along(targets)), collapse = ","), "FROM", source))
  scales = pmax(1, as.numeric(scale_row[1L, ]))
  residual = create(paste("SELECT", paste(c(fe, "w", sprintf("(%s)/(%s) AS %s", targets, dbglm_number(scales), names)), collapse = ","), "FROM", source))
  dictionaries = lapply(fe, function(f) create(paste("SELECT", f, "AS id,", paste0("CAST(0.0 AS DOUBLE) AS ", names, collapse = ","), "FROM", source, "GROUP BY", f)))
  dictionary_spare = rep(list(NULL), 2L)
  means = rep(list(NULL), 2L)
  spare_residual = NULL
  error = Inf
  for (iteration in seq_len(max_iter)) {
    for (i in seq_along(fe)) {
      f = fe[[i]]
      next_mean = create(paste(
        "SELECT", f, "AS id,", paste(paste0(sum_sql(paste0("w*", names)), "/", sum_sql("w"), " AS ", names), collapse = ","),
        "FROM", residual, "GROUP BY", f
      ), target = means[[i]])
      next_dictionary = create(paste(
        "SELECT d.id,", paste(sprintf("d.%s+m.%s AS %s", names, names, names), collapse = ","),
        "FROM", dictionaries[[i]], "d JOIN", next_mean, "m ON d.id=m.id"
      ), target = dictionary_spare[[i]])
      next_residual = create(paste(
        "SELECT", paste(c(
          paste0("s.", fe), "s.w",
          sprintf("s.%s-m.%s AS %s", names, names, names)
        ), collapse = ","), "FROM", residual,
        "s JOIN", next_mean, "m ON", paste0("s.", f, "=m.id")
      ), target = spare_residual)
      old_residual = residual
      residual = next_residual
      spare_residual = old_residual
      old_dictionary = dictionaries[[i]]
      dictionaries[[i]] = next_dictionary
      dictionary_spare[[i]] = old_dictionary
      if (!native) {
        drop(old_residual)
        drop(old_dictionary)
        if (!is.null(means[[i]])) drop(means[[i]])
      }
      means[[i]] = next_mean
    }
    terms = paste(paste0("ABS(", sum_sql(paste0("w*", names)), "/", sum_sql("w"), ") AS ", names), collapse = ",")
    branches = vapply(fe, function(f) paste("SELECT", terms, "FROM", residual, "GROUP BY", f), character(1))
    diagnostic = backend_get_query(conn, paste0(
      "SELECT ", paste(sprintf("MAX(%s) AS %s", names, names), collapse = ","),
      " FROM (", paste(branches, collapse = " UNION ALL "), ") errors"
    ))
    error = max(as.numeric(diagnostic[1L, ]))
    if (is.finite(error) && error <= tol) break
  }
  if (!is.finite(error) || error > tol) stop("Two-FE prediction recovery did not converge.")
  list(
    dictionaries = dictionaries, residual = residual, scales = scales,
    columns = names, iterations = iteration, error = error
  )
}

# Component labels are always IDs from the first FE. Pointer compression
# shortens propagation while every update stays inside its graph component.
prediction_two_fe_components = function(conn, source, fe, create, drop, max_iter = 100L) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  edges = create(paste("SELECT DISTINCT", fe[[1L]], "AS a,", fe[[2L]], "AS b FROM", source))
  left = create(paste("SELECT a AS id,a AS component FROM", edges, "GROUP BY a"))
  right = create(paste("SELECT b AS id,MIN(a) AS component FROM", edges, "GROUP BY b"))
  spare_left = spare_right = spare_compressed = NULL
  for (iteration in seq_len(max_iter)) {
    next_left = create(paste(
      "SELECT e.a AS id,MIN(r.component) AS component FROM", edges,
      "e JOIN", right, "r ON e.b=r.id GROUP BY e.a"
    ), target = spare_left)
    compressed = create(paste("SELECT a.id,b.component FROM", next_left, "a JOIN", next_left, "b ON a.component=b.id"), target = spare_compressed)
    next_right = create(paste(
      "SELECT e.b AS id,MIN(l.component) AS component FROM", edges,
      "e JOIN", compressed, "l ON e.a=l.id GROUP BY e.b"
    ), target = spare_right)
    old_left = left
    old_right = right
    left = compressed
    right = next_right
    spare_left = old_left
    spare_right = old_right
    spare_compressed = next_left
    if (!native) {
      drop(old_left)
      drop(old_right)
      drop(next_left)
    }
    mismatch = backend_get_query(conn, paste(
      "SELECT 1 AS mismatch FROM", edges, "e JOIN", left,
      "l ON e.a=l.id JOIN", right, "r ON e.b=r.id WHERE l.component<>r.component LIMIT 1"
    ))
    if (!nrow(mismatch)) break
  }
  if (nrow(mismatch)) stop("Two-FE component labels did not converge.")
  list(left = left, right = right, iterations = iteration)
}

# A single FE block has a diagonal nuisance system. Keep its dictionary and
# within-group regressor means in SQL; only a p-by-p estimability system moves
# into R. Two FE blocks use SQL projections; varying slopes/three or more
# blocks still require the portable solver.
build_database_prediction_state = function(
  conn, working, retained, fe, components, x,
  design_names, center, create, n, dictionary_components, beta, expected_rss, constant
) {
  backend_info = detect_backend(conn)
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  p = length(x)
  if ((p + 1)^2 > prediction_state_limit()) prediction_state_overflow()
  if (length(beta) != p || any(!is.finite(beta))) stop("Database prediction requires finite fitted slopes.", call. = FALSE)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  metadata = prediction_key_metadata(conn, retained, dictionary_components[[1L]])
  mx = if (p) paste0("m", seq_len(p)) else character()
  slope = if (p) paste(sprintf("(%s)*s.%s", dbglm_number(beta), x), collapse = "+") else "0.0"
  target = "(s.eta-s.__dbglm_o)"
  fields = c(
    paste0("s.", fe, " AS id"),
    paste0(sum_sql(paste0("s.__dbglm_a*(", target, "-(", slope, "))")), "/", sum_sql("s.__dbglm_a"), " AS gamma"),
    if (p) paste0(sum_sql(paste0("s.__dbglm_a*s.", x)), "/", sum_sql("s.__dbglm_a"), " AS ", mx)
  )
  aggregate = paste("SELECT", paste(fields, collapse = ","), "FROM", working, "s GROUP BY", paste0("s.", fe))
  keys = paste0("k", seq_along(dictionary_components[[1L]]))
  dictionary_keys = paste("SELECT DISTINCT", paste(c(
    paste(fe, "AS id"),
    paste(q(dictionary_components[[1L]]), "AS", keys)
  ), collapse = ","), "FROM", retained)
  dictionary = create(paste0(
    "SELECT a.*,", paste0("d.", keys, collapse = ","),
    " FROM (", aggregate, ") a JOIN (", dictionary_keys, ") d ON a.id=d.id"
  ), materialize = TRUE)
  join = paste("FROM", working, "s JOIN", dictionary, "g ON", paste0("s.", fe, "=g.id"))
  residual = paste0("((", slope, ")+g.gamma-", target, ")")
  diagnostics = if (is.null(expected_rss)) paste0(
    "MAX(ABS(", residual,
    ")) AS error,MAX(ABS(", target, ")) AS scale"
  ) else
    paste0(sum_sql(paste0("s.__dbglm_a*POWER(", residual, ",2)")), " AS rss,", sum_sql(paste0("s.__dbglm_a*POWER(", target, ",2)")), " AS scale")
  diagnostic = backend_get_query(conn, paste("SELECT", diagnostics, join))
  recovery_error = if (is.null(expected_rss)) as.numeric(diagnostic[["error"]]) else abs(as.numeric(diagnostic[["rss"]]) - expected_rss)
  tolerance = if (is.null(expected_rss)) 1e-7 else 1e-6
  scale = max(1, as.numeric(diagnostic[["scale"]]), if (!is.null(expected_rss)) abs(expected_rss))
  if (!is.finite(recovery_error) || recovery_error > tolerance * scale) {
    stop("Database prediction recovery does not reproduce the fitted predictor or RSS.", call. = FALSE)
  }
  null = NULL
  if (p) {
    within = paste0("(", paste("SELECT s.__dbglm_a AS w,", paste(sprintf("(s.%s-g.%s) AS %s", x, mx, x), collapse = ","), join), ")")
    H = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", within)), p)
    raw = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ","), "FROM", working)), p)
    null = prediction_within_null(H, raw)
  }
  list(
    version = 1L, specification = list(
      design_names = design_names, design_center = unname(center),
      dictionaries = list(list(
        columns = unique(components[[1L]]), types = metadata[["types"]],
        sql_types = unname(metadata[["sql_types"]]), table = dictionary, keys = keys, means = mx
      ))
    ),
    coefficients = list(slopes = unname(beta)), sample = list(nobs = n), inference = NULL,
    resources = list(tables = gsub("[\"`]", "", dictionary), storage = "database", version = 1L),
    prediction = list(
      storage = "database", null_basis = null, tolerance = 1e-7,
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant
    )
  )
}

finalize_retained_resources = function(result, resources) {
  retained = result[["fitted_state"]][["resources"]]
  if (identical(retained[["storage"]], "database")) {
    result[["resource_state"]] = retain_fit_resources(resources, retained[["tables"]], retained[["views"]])
  }
  result
}

# Match only the requested rows, in bounded batches. The SQL dictionary keeps
# original key types; incoming integers use exact decimal text, never doubles.
prediction_database_lookup = function(state, lease, data) {
  conn = fitted_resource_connection(lease)
  dictionaries = state[["specification"]][["dictionaries"]]
  result = prediction_database_lookup_one(conn, dictionaries[[1L]], data)
  if (length(dictionaries) > 1L) {
    values = c("gamma", dictionaries[[1L]][["means"]])
    for (i in seq.int(2L, length(dictionaries))) {
      next_lookup = prediction_database_lookup_one(conn, dictionaries[[i]], data)
      result[values] = result[values] + next_lookup[values]
      identified = !is.na(result[["component"]]) & !is.na(next_lookup[["component"]]) &
        result[["component"]] == next_lookup[["component"]]
      result[["gamma"]][!identified] = NA_real_
    }
  }
  result
}

prediction_database_lookup_one = function(conn, dictionary, data) {
  keys = dictionary[["keys"]]
  values = Map(function(name, type) {
    encoded = prediction_key_part(data[[name]], type)
    sub("^[0-9]+:", "", encoded)
  }, dictionary[["columns"]], dictionary[["types"]])
  fields = c("gamma", dictionary[["means"]], dictionary[["component"]])
  result = as.data.frame(stats::setNames(rep(list(rep(NA_real_, nrow(data))), length(fields)), fields))
  if (!nrow(data)) return(result)
  batch = max(1L, min(4096L, floor(prediction_state_limit() / (length(fields) + length(keys) + 1L))))
  casts = vapply(seq_along(keys), function(i) sql_prediction_key_cast(
    conn, paste0("v.", keys[[i]]),
    dictionary[["types"]][[i]], dictionary[["sql_types"]][[i]]
  ), character(1))
  join = paste(sprintf("g.%s=%s", keys, casts), collapse = " AND ")
  format = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "json" else "values"
  for (first in seq.int(1L, nrow(data), by = batch)) {
    rows = seq.int(first, min(nrow(data), first + batch - 1L))
    requested = data.frame(i = as.numeric(rows))
    requested[keys] = lapply(values, function(value) value[rows])
    values_sql = sql_values_frame(conn, requested, max_bytes = 2 * 1024^2, format = format)
    if (is.null(values_sql)) {
      stop("Prediction key batch exceeds the SQL transfer limit.", call. = FALSE)
    }
    source = paste0("(", values_sql, ") v")
    sql = paste0(
      "SELECT v.i,", paste0("g.", fields, collapse = ","), " FROM ", source,
      " LEFT JOIN ", dictionary[["table"]], " g ON ", join, " ORDER BY v.i"
    )
    found = backend_get_query(conn, sql)
    if (nrow(found) != length(rows) || !identical(as.numeric(found[["i"]]), as.numeric(rows))) {
      stop("Prediction dictionary no longer has unique fitted keys.", call. = FALSE)
    }
    result[rows, ] = lapply(found[fields], as.numeric)
  }
  result
}

predict_from_database_state = function(state, design, lookup) {
  if (is.null(lookup)) stop("Database prediction requires its fitted key lookup.", call. = FALSE)
  specification = state[["specification"]]
  x = sweep(
    design[, specification[["design_names"]], drop = FALSE], 2L,
    specification[["design_center"]], "-"
  )
  prediction = prediction_matrix_product(x, state[["coefficients"]][["slopes"]]) + lookup[["gamma"]] + state[["prediction"]][["constant"]]
  null = state[["prediction"]][["null_basis"]]
  if (!is.null(null)) {
    means = as.matrix(lookup[specification[["dictionaries"]][[1L]][["means"]]])
    scale = state[["prediction"]][["design_scale"]]
    if (!is.null(scale)) {
      x = sweep(x, 2L, scale, "/")
      means = sweep(means, 2L, scale, "/")
    }
    score = (x - means) %*% null
    norm = sqrt(1 + rowSums(x * x) + rowSums(means * means))
    prediction[sqrt(rowSums(score * score)) > state[["prediction"]][["tolerance"]] * norm] = NA_real_
  }
  prediction[!is.finite(prediction)] = NA_real_
  prediction
}

# Shared bounded rank certificate for one- and two-FE SQL dictionaries.
prediction_within_null = function(H, raw) {
  p = nrow(H)
  if (!p) return(NULL)
  absorbed = diag(H) <= 1e-12 * pmax(diag(raw), .Machine[["double.xmin"]])
  H[absorbed, ] = 0
  H[, absorbed] = 0
  system = list(
    n = p + 1L, big = 1L, rest = 1L + seq_len(p), Dg = 1,
    Cm = matrix(0, 1L, p), E = H, diag = c(1, diag(raw))
  )
  solver = fe_solver(system, new.env(parent = emptyenv()), retain_null = TRUE)
  if (is.null(solver[["null_basis"]])) NULL else solver[["null_basis"]][-1L, , drop = FALSE]
}

build_database_two_fe_state = function(
  conn, working, retained, fe, components, x,
  design_names, center, create, drop, n, dictionary_components, beta, expected_rss,
  constant, recovery_control
) {
  backend_info = detect_backend(conn)
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  p = length(x)
  if ((p + 1)^2 > prediction_state_limit()) prediction_state_overflow()
  if (length(beta) != p || any(!is.finite(beta))) stop("Database prediction requires finite fitted slopes.", call. = FALSE)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  slope = if (p) paste(sprintf("(%s)*s.%s", dbglm_number(beta), x), collapse = "+") else "0.0"
  target = "__dbglm_prediction_target"
  fields = c(paste0("s.", fe), "s.__dbglm_a AS w", paste0("(s.eta-s.__dbglm_o-(", slope, ")) AS ", target), if (p) paste0("s.", x))
  source = create(paste("SELECT", paste(fields, collapse = ","), "FROM", working, "s"))
  projection = prediction_two_fe_projection(conn, source, fe, c(target, x), create, drop,
    tol = min(recovery_control[["tol"]], 1e-10), max_iter = recovery_control[["max_iter"]]
  )
  scale = projection[["scales"]]
  within = projection[["residual"]]
  columns = projection[["columns"]]
  diagnostics = if (is.null(expected_rss)) paste0("MAX(ABS(", columns[[1L]], "))*(", dbglm_number(scale[[1L]]), ") AS error") else
    paste(sum_sql(paste0("w*POWER(", columns[[1L]], "*(", dbglm_number(scale[[1L]]), "),2)")), "AS rss")
  diagnostic = backend_get_query(conn, paste("SELECT", diagnostics, "FROM", within))
  recovery_error = if (is.null(expected_rss)) as.numeric(diagnostic[["error"]]) else abs(as.numeric(diagnostic[["rss"]]) - expected_rss)
  if (is.null(expected_rss)) {
    bound = 1e-7 * scale[[1L]]
  } else {
    norm = backend_get_query(conn, paste("SELECT", sum_sql(paste0("w*POWER( ", target, " ,2)")), "AS norm FROM", source))[["norm"]]
    bound = 1e-6 * max(1, abs(expected_rss), as.numeric(norm))
  }
  if (!is.finite(recovery_error) || recovery_error > bound) {
    stop("Database two-FE recovery does not reproduce the fitted predictor or RSS.", call. = FALSE)
  }
  null = NULL
  if (p) {
    H = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns[-1L], "w"), collapse = ","), "FROM", within)), p)
    normalized = paste0("(SELECT w,", paste(sprintf("(%s)/(%s) AS %s", x, dbglm_number(scale[-1L]), columns[-1L]), collapse = ","), " FROM ", source, ")")
    raw = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns[-1L], "w"), collapse = ","), "FROM", normalized)), p)
    null = prediction_within_null(H, raw)
  }
  labels = prediction_two_fe_components(conn, source, fe, create, drop,
    max_iter = recovery_control[["max_iter"]]
  )
  component_tables = list(labels[["left"]], labels[["right"]])
  mx = if (p) paste0("m", seq_len(p)) else character()
  dictionaries = lapply(seq_along(fe), function(i) {
    metadata = prediction_key_metadata(conn, retained, dictionary_components[[i]])
    keys = paste0("k", seq_along(dictionary_components[[i]]))
    dictionary_keys = paste("SELECT DISTINCT", paste(c(
      paste(fe[[i]], "AS id"),
      paste(q(dictionary_components[[i]]), "AS", keys)
    ), collapse = ","), "FROM", retained)
    fields = c(
      "a.id", paste0("a.", columns[[1L]], "*(", dbglm_number(scale[[1L]]), ") AS gamma"),
      sprintf("a.%s*(%s) AS %s", columns[-1L], dbglm_number(scale[-1L]), mx),
      "l.component", paste0("d.", keys)
    )
    dictionary = create(paste(
      "SELECT", paste(fields, collapse = ","), "FROM", projection[["dictionaries"]][[i]],
      "a JOIN", component_tables[[i]], "l ON a.id=l.id JOIN (", dictionary_keys, ") d ON a.id=d.id"
    ))
    list(
      columns = unique(components[[i]]), types = metadata[["types"]], sql_types = unname(metadata[["sql_types"]]),
      table = dictionary, keys = keys, means = mx, component = "component"
    )
  })
  list(
    version = 1L, specification = list(design_names = design_names, design_center = unname(center), dictionaries = dictionaries),
    coefficients = list(slopes = unname(beta)), sample = list(nobs = n), inference = NULL,
    resources = list(tables = gsub("[\"`]", "", vapply(dictionaries, `[[`, character(1), "table")), storage = "database", version = 1L),
    prediction = list(
      storage = "database", null_basis = null, tolerance = 1e-7, design_scale = scale[-1L],
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant,
      iterations = projection[["iterations"]], projection_error = projection[["error"]], component_iterations = labels[["iterations"]]
    )
  )
}

# Materialized, explicitly keyed output. No input ordering is assumed and no
# prediction panel is implicitly collected into R.
predict_sql = function(object, newdata, conn, table, path, row_id, prediction_table, type = "link") {
  validate_row_id_arg(row_id)
  if (is.null(row_id)) stop("Database prediction requires explicit `row_id` columns.", call. = FALSE)
  if (inherits(object, "dbivreg")) stop("Database prediction output is not yet supported for IV models.", call. = FALSE)
  if (identical(object[["strategy"]], "mundlak")) stop("Database prediction output is not yet supported for Mundlak models.", call. = FALSE)
  if (sum(!vapply(list(newdata, table, path), is.null, logical(1))) != 1L) {
    stop("Provide exactly one of `newdata`, `table`, or `path` for database prediction.", call. = FALSE)
  }
  if (!is.null(prediction_table) && (!is.character(prediction_table) || length(prediction_table) != 1L ||
    is.na(prediction_table) || !nzchar(prediction_table))) stop("`prediction_table` must be a single nonempty table name.", call. = FALSE)
  caller_connection = !is.null(conn)
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) fitted_resource_connection(inherited)
  state = object[["fitted_state"]]
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) {
    owner = fitted_resource_connection(object[["resource_state"]])
    if (!is.null(conn) && !fitted_resource_same_connection(conn, owner)) {
      stop("Database-backed prediction requires its fitted owning session.", call. = FALSE)
    }
    conn = owner
  }
  if (!is.null(prediction_table) && (!caller_connection ||
    (!is.null(inherited) && isTRUE(inherited[["session"]][["own_conn"]])) ||
    (database && isTRUE(object[["resource_state"]][["session"]][["own_conn"]])))) {
    stop("An explicit prediction export requires a caller-owned connection.", call. = FALSE)
  }
  required = prediction_required_columns(object)
  labels = NULL
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) stop("`newdata` must be a data frame.", call. = FALSE)
    incoming = prediction_dataframe_source(object, newdata, unique(c(required, row_id)))
    newdata = incoming[["data"]]
    labels = incoming[["labels"]]
  }
  resources = fit_resources()
  on.exit(release_fit_resources(resources), add = TRUE)
  if (database) resources[["retained_session"]] = object[["resource_state"]][["session"]]
  # Upload incoming R values as typed text on both engines. This also avoids
  # bigint="numeric" registrations reinterpreting integer64 storage bits.
  if (!is.null(newdata)) {
    if (is.null(conn)) {
      conn = DBI::dbConnect(duckdb::duckdb(bigint = "integer64"), shutdown = TRUE)
      resources[["own_conn"]] = TRUE
    }
    resources[["conn"]] = conn
    inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
    if (!is.null(inherited)) fitted_resource_connection(inherited)
    from_statement = paste("FROM", prediction_upload_frame(conn, newdata, resources))
  } else {
    if (is.null(conn) && !is.null(table) && !inherits(table, "tbl_lazy")) stop("A table prediction source requires `conn`.", call. = FALSE)
    setup = setup_db_connection(conn, table, NULL, path, caller = "prediction", resources = resources)
    conn = setup[["conn"]]
    from_statement = setup[["from_statement"]]
  }
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Database prediction output currently supports DuckDB and native Databricks.", call. = FALSE)
  }
  if (!is.null(prediction_table) && sample_table_exists(conn, prediction_table)) stop("Prediction export table already exists.", call. = FALSE)
  context = sql_prediction_context(object, conn, labels)
  source = sql_prediction_source(object, context, conn, from_statement, row_id, required, labels)
  dictionaries = prediction_portable_tables(context, conn, resources)
  prediction = safe_identifiers("prediction", taken = row_id, fold_case = folds_identifier_case(conn))
  query = sql_prediction_equation(object, context, source, dictionaries, conn, row_id, prediction)
  create = function(query) {
    name = basename(tempfile("dbreg_prediction_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, sql_sample_table(conn, name, backend_sql(conn, query), c(row_id, prediction)))
    name
  }
  output = create(query)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  validate_source_row_id(conn, paste("FROM", q(output)), row_id)
  if (inherits(object, "dbglm") && type == "response") {
    output = create(sql_prediction_response(conn, output, row_id, prediction, object[["family"]][["link"]]))
  }
  if (is.null(prediction_table)) {
    lease = retain_fit_resources(resources, output)
  } else {
    DBI::dbExecute(conn, sql_sample_table(conn, prediction_table, paste("SELECT * FROM", q(output)), c(row_id, prediction)))
    output = prediction_table
    lease = retain_fit_resources(resources, character())
  }
  structure(list(
    conn = fitted_resource_public_connection(lease), table = output, row_id = row_id,
    prediction_column = prediction, resource_state = lease
  ), class = "dbreg_prediction")
}

#' Database prediction results
#'
#' Keyed prediction output stays in its database session. `head()` collects a
#' bounded preview; `as.data.frame()` explicitly collects the complete result.
#' `close()` releases owned output tables. Caller-named exports survive release.
#' Output leases share their owning session, so closing a fitted model does not
#' invalidate its prediction results. A model-owned session closes after its
#' final lease is released. Live results cannot be restored by serialization.
#'
#' @param x,con A `dbreg_prediction` result.
#' @param n Number of rows to preview.
#' @param ... Additional arguments are not supported.
#' @return `as.data.frame()` and `head()` return a data frame containing source
#'   keys and the prediction column. `close()` returns `NULL` invisibly.
#' @name dbreg-prediction
NULL

collect_sql_prediction = function(x, n = NULL) {
  conn = fitted_resource_connection(x[["resource_state"]])
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  from = paste("FROM", q(x[["table"]]), if (!is.null(n)) paste("LIMIT", format(n, scientific = FALSE)) else "")
  # Inspect types before collection so exact keys are cast in the first fetch.
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", q(x[["table"]]), "LIMIT 0"))
  keys = x[["row_id"]]
  types = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(keys), ") AS t", seq_along(keys), collapse = ","), "FROM", q(x[["table"]]), "LIMIT 1"))
  if (!nrow(types)) return(schema)
  exact_text = vapply(seq_along(keys), function(i) {
    type = tolower(types[[i]][[1L]])
    if (startsWith(type, "decimal")) return(TRUE)
    value = schema[[keys[[i]]]]
    if (prediction_sql_kind(type) != "integer" || !is.numeric(value) || inherits(value, "integer64")) return(FALSE)
    range = DBI::dbGetQuery(conn, paste0(
      "SELECT MAX(ABS(CAST(", q(keys[[i]]), " AS DOUBLE))) AS maximum FROM (SELECT ",
      q(keys[[i]]), " ", from, ") dbreg_key_range"
    ))[["maximum"]]
    isTRUE(as.numeric(range) >= 2^53)
  }, logical(1))
  fields = names(schema)
  select = q(fields)
  if (any(exact_text)) {
    text = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "STRING" else "VARCHAR"
    index = match(keys[exact_text], fields)
    select[index] = paste0("CAST(", q(fields[index]), " AS ", text, ") AS ", q(fields[index]))
  }
  output = DBI::dbGetQuery(conn, paste("SELECT", paste(select, collapse = ","), from))
  output
}

#' @rdname dbreg-prediction
#' @export
as.data.frame.dbreg_prediction = function(x, ...) {
  check_unused_args(list(...))
  collect_sql_prediction(x)
}

#' @rdname dbreg-prediction
#' @importFrom utils head
#' @export
head.dbreg_prediction = function(x, n = 6L, ...) {
  check_unused_args(list(...))
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) stop("`n` must be a nonnegative whole number.", call. = FALSE)
  collect_sql_prediction(x, n)
}

#' @rdname dbreg-prediction
#' @export
print.dbreg_prediction = function(x, ...) {
  check_unused_args(list(...))
  cat("Database predictions in ", x[["table"]], "; keys: ", paste(x[["row_id"]], collapse = ", "),
    "; prediction column: ", x[["prediction_column"]], ".\n",
    sep = ""
  )
  invisible(x)
}

#' @rdname dbreg-prediction
#' @export
close.dbreg_prediction = function(con, ...) {
  check_unused_args(list(...))
  fitted_resource_release(con[["resource_state"]])
}

prediction_i_specs = function(object) {
  if (inherits(object, "dbglm")) object[["model_spec"]][["i_specs"]] else object[["i_specs"]]
}

prediction_required_columns = function(object) {
  if (inherits(object, "dbglm")) {
    variables = object[["model_spec"]][["source_vars"]]
    for (name in intersect(variables, names(object[["factor_sources"]]))) variables[variables == name] = object[["factor_sources"]][[name]]
  } else variables = all.vars(stats::formula(object[["fml"]], lhs = 0, rhs = 1))
  state = object[["fitted_state"]]
  unique(c(
    variables, object[["offset"]],
    unlist(lapply(state[["specification"]][["dictionaries"]], `[[`, "columns")),
    unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable"))
  ))
}

prediction_factor_sources = function(object) {
  if (inherits(object, "dbglm")) {
    variables = names(object[["model_spec"]][["factor_levels"]])
    for (name in intersect(variables, names(object[["factor_sources"]]))) variables[variables == name] = object[["factor_sources"]][[name]]
    unique(variables)
  } else {
    rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    unique(c(
      intersect(names(object[["factor_levels"]]), all.vars(rhs)),
      unname(explicit_factor_sources(attr(stats::terms(rhs), "term.labels")))
    ))
  }
}

prediction_dataframe_source = function(object, data, required) {
  data = as.data.frame(data)
  absent = setdiff(required, names(data))
  if (length(absent)) stop("newdata is missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  state = object[["fitted_state"]]
  character_keys = unlist(lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
    dictionary[["columns"]][dictionary[["types"]] == "character"]
  }), use.names = FALSE)
  factors = unique(c(prediction_factor_sources(object), character_keys))
  labels = safe_identifiers(paste0("dbreg_prediction_factor", seq_along(factors), "__"), taken = names(data))
  if (!length(factors)) labels = character()
  for (i in seq_along(factors)) data[[labels[[i]]]] = as.character(data[[factors[[i]]]])
  names(labels) = factors
  if (length(prediction_i_specs(object))) {
    eligible = rep(TRUE, nrow(data))
    for (name in prediction_required_columns(object)) {
      value = data[[name]]
      eligible = eligible & !is.na(value)
      if (is.numeric(value) && !inherits(value, "integer64")) eligible = eligible & is.finite(value)
    }
  }
  for (spec in prediction_i_specs(object)) {
    first = spec[["first"]]
    public = if (inherits(object, "dbglm")) object[["factor_sources"]][[first]] else first
    if (is.null(public)) public = first
    index = rep(NA_integer_, nrow(data))
    index[eligible] = i_match_values(spec, data[[public]][eligible], strict = FALSE)
    label = safe_identifiers(paste0("dbreg_prediction_i", length(labels) + 1L, "__"), taken = names(data))
    data[[label]] = spec[["labels"]][index]
    # Leave rejection to the SQL validity/known-FE checks, in the same order
    # as R prediction. The marker cannot alias an actual fitted label.
    marker = "dbreg_prediction_unseen_i__"
    while (marker %in% spec[["labels"]]) marker = paste0(marker, "_")
    data[[label]][eligible & is.na(index)] = marker
    labels[first] = label
  }
  list(data = data[unique(c(required, labels))], labels = labels)
}

sql_prediction_context = function(object, conn, dataframe_labels = NULL) {
  state = object[["fitted_state"]]
  if (is.null(state[["prediction"]])) state = NULL
  if (length(object[["fe"]]) && is.null(state)) {
    stop("Database prediction with fixed effects requires fitting with retain='prediction'.", call. = FALSE)
  }
  if (!is.null(state) && !identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  design = sql_fitted_prediction_design(object, conn, factor_labels = !is.null(dataframe_labels))
  for (name in names(design[["factor_levels"]])) {
    public = design[["sources"]][[name]]
    if (!is.null(dataframe_labels)) design[["sources"]][name] = dataframe_labels[[if (name %in% names(dataframe_labels)) name else public]]
  }
  if (is.null(state)) {
    beta = stats::coef(object)
    if (inherits(object, "dbglm")) names(beta) = object[["design_names"]]
    index = if (!is.null(object[["linear_design"]])) {
      match(object[["linear_design"]][["fitted_ids"]], design[["ids"]])
    } else match(names(beta), design[["names"]])
    if (anyNA(index)) stop("Could not align fitted SQL prediction coefficients.", call. = FALSE)
    design[["expressions"]] = design[["expressions"]][index]
    design[["names"]] = design[["names"]][index]
    return(list(
      kind = "none", design = design, beta = unname(beta), center = rep(0, length(beta)),
      constant = 0, state = NULL, offset = object[["offset"]], varying = character()
    ))
  }
  expected = state[["specification"]][["design_names"]]
  index = match(expected, design[["names"]])
  if (anyNA(index)) stop("Could not align the fitted learned-effect SQL design.", call. = FALSE)
  design[["expressions"]] = design[["expressions"]][index]
  design[["names"]] = expected
  database = identical(state[["prediction"]][["storage"]], "database")
  beta = if (database) state[["coefficients"]][["slopes"]] else
    state[["coefficients"]][["joint"]][state[["coefficients"]][["fe_n"]] + seq_along(expected)]
  constant = state[["prediction"]][["constant"]]
  if (is.null(constant)) constant = 0
  list(
    kind = if (database) "database" else "portable", design = design, beta = beta,
    center = state[["specification"]][["design_center"]], constant = constant,
    state = state, offset = object[["offset"]], varying = unique(unlist(lapply(state[["specification"]][["operations"]], `[[`, "variable")))
  )
}

prediction_portable_tables = function(context, conn, resources) {
  state = context[["state"]]
  if (context[["kind"]] == "none") return(list())
  if (context[["kind"]] == "database") return(state[["specification"]][["dictionaries"]])
  dictionaries = state[["specification"]][["dictionaries"]]
  coefficients = state[["coefficients"]]
  operations = state[["specification"]][["operations"]]
  null = state[["prediction"]][["null_basis"]]
  lapply(seq_along(dictionaries), function(i) {
    dictionary = dictionaries[[i]]
    frame = prediction_decode_keys(dictionary)
    key_types = vapply(seq_along(dictionary[["types"]]), function(j) {
      switch(dictionary[["types"]][[j]],
        integer = prediction_integer_sql_type(frame[[j]], conn),
        numeric = "DOUBLE",
        float = "REAL",
        character = "STRING",
        logical = "BOOLEAN",
        date = "DATE"
      )
    }, character(1))
    types = stats::setNames(key_types, names(frame))
    frame[["id"]] = as.numeric(dictionary[["id"]])
    types[["id"]] = "DOUBLE"
    members = which(vapply(operations, `[[`, integer(1), "group") == i)
    for (op in members) {
      indices = coefficients[["starts"]][[op]] + match(dictionary[["id"]], coefficients[["ids"]][[op]])
      frame[[paste0("gamma", op)]] = coefficients[["joint"]][indices]
      types[[paste0("gamma", op)]] = "DOUBLE"
      if (!is.null(null)) for (j in seq_len(ncol(null))) {
        name = paste0("null", op, "_", j)
        frame[[name]] = null[indices, j]
        types[[name]] = "DOUBLE"
      }
    }
    dictionary[["table"]] = prediction_upload_frame(conn, frame, resources, types)
    dictionary[["keys"]] = paste0("k", seq_along(key_types))
    dictionary[["sql_types"]] = key_types
    dictionary[["operations"]] = members
    dictionary
  })
}

sql_prediction_numeric = function(conn, column, class) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  if (class == "Date") {
    return(if (native) paste0("DATEDIFF(", column, ",DATE '1970-01-01')") else
      paste0("DATE_DIFF('day',DATE '1970-01-01',", column, ")"))
  }
  if (class == "POSIXct") {
    return(if (native) paste0("unix_micros(", column, ")/1e6") else paste0("epoch(", column, ")"))
  }
  if (!class %in% c("numeric", "integer", "integer64", "logical")) stop("Numeric prediction inputs require numeric or date/time SQL columns.", call. = FALSE)
  paste0("CAST(", column, " AS DOUBLE)")
}

sql_prediction_key_cast = function(conn, value, type, sql_type, input_type = "STRING") {
  cast = paste0("TRY_CAST(", value, " AS ", sql_type, ")")
  if (type == "integer" && grepl("^(varchar|char|string|enum)", tolower(input_type))) {
    valid = if (isTRUE(detect_backend(conn)[["is_databricks"]])) paste0(value, " RLIKE '^[+-]?[0-9]+$'") else
      paste0("regexp_full_match(", value, ",'[+-]?[0-9]+')")
    return(paste0("CASE WHEN ", valid, " THEN ", cast, " ELSE CAST(NULL AS ", sql_type, ") END"))
  }
  if (type == "date" && grepl("^(varchar|char|string|enum)", tolower(input_type))) {
    # Date strings match the fitted textual date, rather than silently accepting
    # a timestamp or a differently formatted string for the same calendar day.
    text = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "STRING" else "VARCHAR"
    valid = paste0(
      "((CAST(", cast, " AS ", text, ")=", value, " AND ", value,
      " NOT IN ('infinity','-infinity')) OR ", value, " IN ('Inf','-Inf'))"
    )
    return(paste0("CASE WHEN ", valid, " THEN ", cast, " ELSE CAST(NULL AS ", sql_type, ") END"))
  }
  cast
}

# Recompile from fitted coding, never category discovery in prediction data.
# The same term expanders serve estimation and this SQL prediction path.
sql_fitted_prediction_design = function(object, conn, factor_labels = FALSE) {
  glm = inherits(object, "dbglm")
  state = object[["fitted_state"]]
  if (glm) {
    specification = object[["model_spec"]]
    rhs = specification[["rhs"]]
    variables = specification[["source_vars"]]
    if (any(!is_safe_identifier(variables))) {
      stop("SQL prediction requires SQL-safe fitted GLM regressor names.", call. = FALSE)
    }
    aliases = variables
    sources = stats::setNames(variables, aliases)
    for (name in intersect(names(object[["factor_sources"]]), aliases)) sources[name] = object[["factor_sources"]][[name]]
    levels = specification[["factor_levels"]]
    classes = specification[["factor_classes"]]
    if (factor_labels) classes = stats::setNames(rep("character", length(levels)), names(levels))
    info = list(
      i_specs = specification[["i_specs"]], types = as.list(stats::setNames(ifelse(aliases %in% names(levels), "factor", "numeric"), aliases)),
      levels = levels, factor_classes = classes, level_sql = if (factor_labels)
        lapply(levels, function(value) sql_value_literals(conn, as.character(value))) else sql_factor_literals(conn, levels, classes)
    )
    if (!factor_labels) for (name in names(specification[["i_specs"]])) {
      values = stats::setNames(list(specification[["i_specs"]][[name]][["values"]]), name)
      info[["level_sql"]][[name]] = sql_factor_literals(conn, values, classes)[[name]]
    }
    info[["classes"]] = specification[["source_classes"]]
    info[["numeric_sql"]] = stats::setNames(sql_cast_double(aliases, conn), aliases)
    design = sql_model_matrix_from_info(attr(stats::terms(rhs), "term.labels"), info,
      intercept = attr(stats::terms(rhs), "intercept") == 1L || length(object[["fe"]]) > 0L
    )
    names = design[["coef_names"]]
    expressions = design[["select_exprs"]]
    if (!identical(names, specification[["columns"]][["name"]])) {
      stop("Could not reconstruct the fitted GLM SQL design.", call. = FALSE)
    }
    intercept = isTRUE(object[["has_intercept"]])
  } else {
    rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    variables = all.vars(rhs)
    aliases = safe_identifiers(variables, fold_case = folds_identifier_case(conn))
    sources = stats::setNames(variables, aliases)
    renamed = aliases != variables
    source_map = stats::setNames(variables[renamed], aliases[renamed])
    rhs = rename_expr_vars(rhs, stats::setNames(aliases, variables))
    source_aliases = stats::setNames(aliases, variables)
    raw_rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    forced = explicit_factor_sources(attr(stats::terms(raw_rhs), "term.labels"))
    relevant = c(variables, names(forced))
    levels = object[["factor_levels"]][intersect(names(object[["factor_levels"]]), relevant)]
    if (length(levels)) names(levels) = restore_term_names(names(levels), source_aliases)
    classes = object[["factor_classes"]][intersect(names(object[["factor_classes"]]), relevant)]
    if (length(classes)) names(classes) = restore_term_names(names(classes), source_aliases)
    forced_safe = explicit_factor_sources(attr(stats::terms(rhs), "term.labels"))
    for (term in names(forced_safe)) sources[[term]] = variables[match(forced_safe[[term]], aliases)]
    if (length(forced)) {
      unused = setdiff(unname(forced), model_vars(raw_rhs, skip_factors = TRUE))
      unused = aliases[match(unused, variables)]
      sources = sources[!names(sources) %in% unused]
      levels = levels[!names(levels) %in% unused]
      classes = classes[!names(classes) %in% unused]
    }
    if (factor_labels) classes = stats::setNames(rep("character", length(levels)), names(levels))
    # Scalar FE terms affect the native interaction reference rule; composite
    # FE ranks are distinct from their component variables in that rule.
    scalar_fe = unlist(lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
      if (length(dictionary[["columns"]]) == 1L) dictionary[["columns"]] else character()
    }), use.names = FALSE)
    fe = aliases[match(intersect(scalar_fe, variables), variables)]
    info = list(
      types = as.list(stats::setNames(ifelse(aliases %in% names(levels), "factor", "numeric"), aliases)),
      levels = levels, factor_classes = classes, level_sql = if (factor_labels)
        lapply(levels, function(value) sql_value_literals(conn, as.character(value))) else sql_factor_literals(conn, levels, classes)
    )
    info[["i_specs"]] = lapply(object[["i_specs"]], function(spec) {
      spec[["first"]] = source_aliases[[spec[["first"]]]]
      if (!is.null(spec[["second"]])) spec[["second"]] = source_aliases[[spec[["second"]]]]
      spec[["term"]] = restore_term_names(spec[["term"]], source_aliases)
      spec
    })
    for (spec in info[["i_specs"]]) {
      first = spec[["first"]]
      if (!is.null(spec[["source_type"]])) info[["types"]][[first]] = spec[["source_type"]]
      info[["levels"]][[first]] = if (factor_labels) spec[["labels"]] else spec[["values"]]
      info[["level_sql"]][[first]] = if (factor_labels) sql_value_literals(conn, spec[["labels"]]) else
        sql_factor_literals(conn, stats::setNames(list(spec[["values"]]), first), classes)[[first]]
      levels[[first]] = if (factor_labels) spec[["labels"]] else spec[["values"]]
    }
    info[["classes"]] = stats::setNames(rep("numeric", length(aliases)), aliases)
    info[["classes"]][names(classes)] = classes
    for (spec in info[["i_specs"]]) {
      second = spec[["second"]]
      if (is.null(second)) next
      info[["classes"]][second] = spec[["second_class"]]
      if (identical(spec[["second_class"]], "logical") && factor_labels) {
        info[["logical_sql"]][[second]] = paste0(second, " = ", sql_value_literals(conn, "TRUE"))
      }
    }
    info[["factor_sources"]] = as.list(stats::setNames(as.character(DBI::dbQuoteIdentifier(conn, as.character(names(forced_safe)))), names(forced_safe)))
    info[["numeric_sql"]] = stats::setNames(sql_cast_double(aliases, conn), aliases)
    intercept = attr(stats::terms(rhs), "intercept") == 1L
    design = sql_model_matrix_from_info(attr(stats::terms(rhs), "term.labels"), info,
      fe_vars = fe, intercept = intercept, fold_case = folds_identifier_case(conn)
    )
    inputs = new.env(parent = emptyenv())
    inputs[["source_map"]] = source_map
    inputs[["name_map"]] = source_map
    register_factor_display(inputs, design[["factor_levels"]], design[["level_aliases"]],
      design[["col_names"]], design[["aliases"]],
      i_columns = design[["i_columns"]]
    )
    map = coefficient_name_map(inputs[["name_map"]], source_map)
    names = standardize_coef_names(restore_names(design[["aliases"]], map))
    expressions = design[["select_exprs"]]
    intercept = intercept && !length(object[["fe"]])
  }
  if (intercept) {
    expressions = c("1.0", expressions)
    names = c("(Intercept)", names)
  }
  list(
    expressions = expressions, names = names, sources = sources,
    ids = if (!glm) c(if (intercept) "intercept", design[["columns"]][["id"]]),
    factor_levels = levels, factor_classes = classes,
    level_sql = info[["level_sql"]],
    factor_labels = factor_labels,
    numeric_sources = setdiff(aliases, names(levels))
  )
}

prediction_sql_kind = function(type) {
  type = tolower(gsub(" ", "", type, fixed = TRUE))
  if (grepl("^(u?(tinyint|smallint|int|integer|bigint|hugeint)|decimal\\([0-9]+,0\\))$", type)) "integer" else if (grepl("^(double|float|real|decimal)", type)) "numeric" else if (grepl("^(varchar|char|string|enum)", type)) "character" else if (type %in% c("boolean", "bool")) "logical" else if (type == "date") "date" else if (startsWith(type, "timestamp")) "datetime" else "unsupported"
}

prediction_factor_class = function(design, name) {
  class = unname(design[["factor_classes"]][name])
  if (!length(class) || is.na(class)) class = class(design[["factor_levels"]][[name]])[[1L]]
  class
}

sql_prediction_source = function(object, context, conn, from_statement, row_id, raw_required, dataframe_labels = NULL) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  dictionaries = context[["state"]][["specification"]][["dictionaries"]]
  key_sources = lapply(dictionaries, function(dictionary) {
    columns = dictionary[["columns"]]
    for (j in seq_along(columns)) {
      if (dictionary[["types"]][[j]] == "character" && !is.null(dataframe_labels)) columns[[j]] = dataframe_labels[[columns[[j]]]]
    }
    columns
  })
  needed = unique(c(raw_required, unname(context[["design"]][["sources"]]), unlist(key_sources), row_id))
  schema = DBI::dbGetQuery(conn, paste("SELECT *", from_statement, "LIMIT 0"))
  absent = setdiff(needed, names(schema))
  if (length(absent)) stop("Prediction source is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  type_row = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(needed), ") AS t", seq_along(needed), collapse = ","), from_statement, "LIMIT 1"))
  types = if (nrow(type_row)) stats::setNames(vapply(type_row, function(value) value[[1L]], character(1)), needed) else
    prediction_frame_types(schema[needed], conn)
  kinds = vapply(types, prediction_sql_kind, character(1))
  classes = vapply(schema[needed], function(value) class(value)[[1L]], character(1))
  reserved = c(names(context[["design"]][["sources"]]), "__dbglm_y", "__dbglm_a", "__dbglm_o", "__prediction_eta")
  raw_aliases = safe_identifiers(paste0("dbreg_prediction_raw", seq_along(needed), "__"),
    taken = reserved,
    fold_case = folds_identifier_case(conn)
  )
  names(raw_aliases) = needed
  key_aliases = safe_identifiers(paste0("dbreg_prediction_id", seq_along(row_id), "__"), taken = c(reserved, raw_aliases), fold_case = folds_identifier_case(conn))
  x = if (length(context[["beta"]])) safe_identifiers(paste0("dbreg_prediction_x", seq_along(context[["beta"]]), "__"),
    taken = c(reserved, raw_aliases, key_aliases), fold_case = folds_identifier_case(conn)
  ) else character()
  varying = context[["varying"]]
  varying_aliases = if (length(varying)) safe_identifiers(paste0("dbreg_prediction_v", seq_along(varying), "__"),
    taken = c(reserved, raw_aliases, key_aliases, x), fold_case = folds_identifier_case(conn)
  ) else character()
  names(varying_aliases) = varying
  offset_alias = safe_identifiers("dbreg_prediction_offset__", taken = c(reserved, raw_aliases, key_aliases, x, varying_aliases), fold_case = folds_identifier_case(conn))
  raw = paste("SELECT", paste(q(needed), "AS", raw_aliases, collapse = ","), from_statement)
  precision = character()
  source_terms = vapply(names(context[["design"]][["sources"]]), function(name) {
    source = context[["design"]][["sources"]][[name]]
    value = paste0("raw.", raw_aliases[[source]])
    if (name %in% context[["design"]][["numeric_sources"]]) {
      value = sql_prediction_numeric(conn, value, classes[[source]])
    } else {
      class = prediction_factor_class(context[["design"]], name)
      kind = kinds[[source]]
      compatible = switch(class,
        character = kind == "character",
        factor = kind == "character",
        ordered = kind == "character",
        logical = kind == "logical",
        Date = kind == "date",
        integer64 = kind %in% c("integer", "numeric"),
        integer = kind %in% c("integer", "numeric"),
        numeric = kind %in% c("integer", "numeric"),
        FALSE
      )
      if (!compatible) stop("Prediction factor `", source, "` requires its fitted column type.", call. = FALSE)
      if (class == "integer64") value = paste0("TRY_CAST(", value, " AS BIGINT)") else if (class %in% c("numeric", "integer")) value = paste0("CAST(", value, " AS DOUBLE)")
    }
    paste(value, "AS", q(name))
  }, character(1))
  # Precision checks use original fields, before a categorical numeric cast.
  for (name in names(context[["design"]][["factor_levels"]])) {
    source = context[["design"]][["sources"]][[name]]
    class = prediction_factor_class(context[["design"]], name)
    value = paste0("r.", raw_aliases[[source]])
    if (class == "integer64" && kinds[[source]] == "numeric") precision = c(
      precision,
      paste0("(ABS(", value, ")>9007199254740992 OR ", value, "<>FLOOR(", value, "))")
    )
    if (class %in% c("numeric", "integer") && kinds[[source]] == "integer") precision = c(
      precision,
      paste0("(", value, ">9007199254740992 OR ", value, "< -9007199254740992)")
    )
  }
  for (i in seq_along(dictionaries)) for (j in seq_along(key_sources[[i]])) {
    source = key_sources[[i]][[j]]
    type = dictionaries[[i]][["types"]][[j]]
    kind = kinds[[source]]
    compatible = switch(type,
      integer = kind %in% c("integer", "numeric", "character"),
      numeric = kind %in% c("integer", "numeric"),
      float = kind %in% c("integer", "numeric"),
      character = kind == "character",
      logical = kind == "logical",
      date = kind %in% c("date", "character"),
      FALSE
    )
    if (!compatible) stop("Prediction FE key `", source, "` requires its fitted column type.", call. = FALSE)
    value = paste0("r.", raw_aliases[[source]])
    if (type == "integer" && kind == "numeric") precision = c(
      precision,
      paste0("(ABS(", value, ")>9007199254740992 OR ", value, "<>FLOOR(", value, "))")
    )
    if (type %in% c("numeric", "float") && kind == "integer") precision = c(
      precision,
      paste0("(", value, ">9007199254740992 OR ", value, "< -9007199254740992)")
    )
  }
  output_keys = vapply(seq_along(row_id), function(i) {
    value = paste0("raw.", raw_aliases[[row_id[[i]]]])
    # SQL NaN is non-NULL on DuckDB/Spark, but is a missing key in R.
    if (kinds[[row_id[[i]]]] == "numeric") value = paste0(
      "CASE WHEN isnan(CAST(", value,
      " AS DOUBLE)) THEN NULL ELSE ", value, " END"
    )
    paste(value, "AS", key_aliases[[i]])
  }, character(1))
  extras = c(
    output_keys,
    vapply(varying, function(name) paste(sql_prediction_numeric(conn, paste0("raw.", raw_aliases[[name]]), classes[[name]]), "AS", varying_aliases[[name]]), character(1)),
    paste(if (is.null(context[["offset"]])) "CAST(0.0 AS DOUBLE)" else
      sql_prediction_numeric(conn, paste0("raw.", raw_aliases[[context[["offset"]]]]), classes[[context[["offset"]]]]), "AS", offset_alias)
  )
  prepared = paste0("(SELECT raw.*,", paste(c(source_terms, extras), collapse = ","), " FROM (", raw, ") raw) r")
  fields = c("r.*", sprintf("(CAST((%s) AS DOUBLE)-(%s)) AS %s", context[["design"]][["expressions"]], dbglm_number(context[["center"]]), x))
  frame = paste0("(SELECT ", paste(fields, collapse = ","), " FROM ", prepared, ") r")
  raw_values = if (length(raw_required)) paste0("r.", raw_aliases[raw_required]) else character()
  valid = if (length(raw_values)) paste0(raw_values, " IS NOT NULL") else character()
  numeric = raw_required[classes[raw_required] %in% c("numeric", "integer", "integer64")]
  if (length(numeric)) {
    values = paste0("r.", raw_aliases[numeric])
    valid = c(valid, sql_is_finite(conn, values))
  }
  list(
    frame = frame, x = x, varying = varying_aliases, offset = offset_alias,
    key_aliases = key_aliases, raw_aliases = raw_aliases, key_sources = key_sources,
    types = types, classes = classes, valid = if (length(valid)) paste(valid, collapse = " AND ") else "TRUE",
    precision = if (length(precision)) paste(precision, collapse = " OR ") else NULL
  )
}

sql_prediction_sum = function(terms) {
  if (length(terms)) paste0("(", paste(terms, collapse = "+"), ")") else "0.0"
}

sql_prediction_equation = function(object, context, source, dictionaries, conn, row_id, prediction) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  product = function(coef, value) sprintf("(%s)*(%s)", dbglm_number(coef), value)
  square = function(value) if (length(value)) paste0("POWER((", value, "),2)") else character()
  x = if (length(source[["x"]])) paste0("r.", source[["x"]]) else character()
  joins = vapply(seq_along(dictionaries), function(i) {
    dictionary = dictionaries[[i]]
    keys = source[["key_sources"]][[i]]
    comparisons = vapply(seq_along(keys), function(j) {
      value = paste0("r.", source[["raw_aliases"]][[keys[[j]]]])
      cast = sql_prediction_key_cast(conn, value, dictionary[["types"]][[j]], dictionary[["sql_types"]][[j]], source[["types"]][[keys[[j]]]])
      paste0("g", i, ".", dictionary[["keys"]][[j]], "=", cast)
    }, character(1))
    paste0(" LEFT JOIN ", dictionary[["table"]], " g", i, " ON ", paste(comparisons, collapse = " AND "))
  }, character(1))
  known = if (length(dictionaries)) paste0("g", seq_along(dictionaries), ".id IS NOT NULL") else "TRUE"
  state = context[["state"]]
  null = state[["prediction"]][["null_basis"]]
  contribution = score = character()
  norm = "1.0"
  if (context[["kind"]] == "database") {
    contribution = paste0("g", seq_along(dictionaries), ".gamma")
    if (length(dictionaries) == 2L) known = c(known, "g1.component=g2.component")
    if (!is.null(null) && ncol(null)) {
      scale = state[["prediction"]][["design_scale"]]
      if (is.null(scale)) scale = rep(1, length(x))
      means = vapply(seq_along(x), function(j) sql_prediction_sum(vapply(seq_along(dictionaries), function(i)
        paste0("g", i, ".", dictionaries[[i]][["means"]][[j]]), character(1))), character(1))
      difference = sprintf("((%s-(%s))/(%s))", x, means, dbglm_number(scale))
      score = vapply(seq_len(ncol(null)), function(j) sql_prediction_sum(product(null[, j], difference)), character(1))
      norm = sql_prediction_sum(c(
        "1.0", square(sprintf("(%s)/(%s)", x, dbglm_number(scale))),
        square(sprintf("(%s)/(%s)", means, dbglm_number(scale)))
      ))
    }
  } else if (context[["kind"]] == "portable") {
    operations = state[["specification"]][["operations"]]
    multipliers = vapply(operations, function(op) if (is.null(op[["variable"]])) "1.0" else
      paste0(
        "(r.", source[["varying"]][[op[["variable"]]]], "-(",
        dbglm_number(if (is.null(op[["center"]])) 0 else op[["center"]]), "))"
      ), character(1))
    contribution = vapply(seq_along(operations), function(j) paste0("(", multipliers[[j]], ")*g", operations[[j]][["group"]], ".gamma", j), character(1))
    if (!is.null(null) && ncol(null)) {
      fe_n = state[["coefficients"]][["fe_n"]]
      score = vapply(seq_len(ncol(null)), function(k) {
        effects = vapply(seq_along(operations), function(j) paste0("(", multipliers[[j]], ")*g", operations[[j]][["group"]], ".null", j, "_", k), character(1))
        sql_prediction_sum(c(product(null[fe_n + seq_along(x), k], x), effects))
      }, character(1))
      norm = sql_prediction_sum(c(square(x), square(multipliers)))
    }
  }
  eta = sql_prediction_sum(c(product(context[["beta"]], x), contribution, dbglm_number(context[["constant"]]), paste0("r.", source[["offset"]])))
  estimable = if (length(score)) paste0(
    sql_prediction_sum(square(score)), "<=",
    dbglm_number(state[["prediction"]][["tolerance"]]^2), "*(", norm, ")"
  ) else "TRUE"
  factor_valid = vapply(names(context[["design"]][["factor_levels"]]), function(name) {
    levels = context[["design"]][["factor_levels"]][[name]]
    class = prediction_factor_class(context[["design"]], name)
    literals = context[["design"]][["level_sql"]][[name]]
    if (!length(literals)) return("FALSE")
    paste0("r.", q(name), " IN (", paste(literals, collapse = ","), ")")
  }, character(1))
  error = function(message) paste0(
    "CAST(", if (isTRUE(detect_backend(conn)[["is_databricks"]])) "raise_error" else "error",
    "(", sql_value_literals(conn, message), ") AS DOUBLE)"
  )
  cases = c(
    paste0("WHEN NOT (", source[["valid"]], ") THEN CAST(NULL AS DOUBLE)"),
    if (!is.null(source[["precision"]])) paste0("WHEN ", source[["precision"]], " THEN ", error("Prediction keys or factor values cannot be represented exactly; supply exact integer columns or integer FE strings.")),
    paste0("WHEN NOT (", paste(known, collapse = " AND "), ") THEN CAST(NULL AS DOUBLE)"),
    if (length(factor_valid)) paste0("WHEN NOT (", paste(factor_valid, collapse = " AND "), ") THEN ", error("newdata contains unseen factor levels.")),
    paste0("WHEN NOT (", estimable, ") THEN CAST(NULL AS DOUBLE)"),
    if (context[["kind"]] != "none") paste0("WHEN NOT ", sql_is_finite(conn, eta), " THEN CAST(NULL AS DOUBLE)")
  )
  output = c(
    paste0("r.", source[["key_aliases"]], " AS ", q(row_id)),
    paste0("CASE ", paste(cases, collapse = " "), " ELSE ", eta, " END AS ", q(prediction))
  )
  paste0("SELECT ", paste(output, collapse = ","), " FROM ", source[["frame"]], paste(joins, collapse = ""))
}

sql_prediction_response = function(conn, table, row_id, prediction, link) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  if (link == "log") {
    # R's Poisson linkinv floors tiny means at epsilon and remains finite above
    # eta=700. The conservative fitting guard must not limit new predictions.
    eta = q(prediction)
    mu = paste0("GREATEST(EXP(", eta, "),", dbglm_number(.Machine[["double.eps"]]), ")")
    return(paste0(
      "SELECT ", paste(q(row_id), collapse = ","), ",CASE WHEN ", eta,
      " IS NULL OR isnan(", eta, ") THEN CAST(NULL AS DOUBLE) ELSE ", mu, " END AS ",
      q(prediction), " FROM ", q(table)
    ))
  }
  keys = paste0("dbreg_prediction_key", seq_along(row_id), "__")
  mean_source = paste0(
    "(SELECT ", paste(q(row_id), "AS", keys, collapse = ","), ",", q(prediction),
    " AS __prediction_eta,CAST(0.0 AS DOUBLE) AS __dbglm_y,CAST(1.0 AS DOUBLE) AS __dbglm_a,",
    "CAST(0.0 AS DOUBLE) AS __dbglm_o FROM ", q(table), ")"
  )
  working = dbglm_working_sql(mean_source, character(), numeric(), link, eta = "__prediction_eta")
  paste0(
    working, "SELECT ", paste(keys, "AS", q(row_id), collapse = ","),
    ",CASE WHEN __prediction_eta IS NULL OR isnan(__prediction_eta) THEN CAST(NULL AS DOUBLE) ELSE mu END AS ",
    q(prediction), " FROM means"
  )
}

prediction_frame_types = function(frame, conn) {
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  vapply(frame, function(value) {
    if (inherits(value, "integer64")) "BIGINT" else if (inherits(value, "Date")) "DATE" else if (inherits(value, "POSIXct")) if (native) "TIMESTAMP" else "TIMESTAMPTZ" else if (is.factor(value) || is.character(value)) if (native) "STRING" else "VARCHAR" else if (is.logical(value)) "BOOLEAN" else if (is.integer(value)) "INTEGER" else if (is.numeric(value) && !is.complex(value)) "DOUBLE" else
      stop("Unsupported prediction input column type: ", class(value)[[1L]], ".", call. = FALSE)
  }, character(1))
}

prediction_frame_text = function(value) {
  if (inherits(value, "integer64")) return(as.character(value))
  if (inherits(value, "POSIXct")) return(format(value, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  if (is.numeric(value) && !inherits(value, "Date")) {
    text = dbglm_number(value)
    text[is.na(value)] = NA_character_
    text[!is.na(value) & value == Inf] = "Infinity"
    text[!is.na(value) & value == -Inf] = "-Infinity"
    return(text)
  }
  as.character(value)
}

# Bounded incoming batches/portable dictionaries, never a collected SQL panel.
prediction_upload_frame = function(conn, frame, resources, sql_types = NULL) {
  if (is.null(sql_types)) sql_types = prediction_frame_types(frame, conn)
  stopifnot(length(sql_types) == ncol(frame), ncol(frame) > 0L)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  name = basename(tempfile("dbreg_prediction_input_"))
  resources[["tables"]] = c(resources[["tables"]], name)
  fields = sprintf("CAST(%s AS %s) AS %s", q(names(frame)), sql_types, q(names(frame)))
  empty = paste("SELECT", paste(sprintf("CAST(NULL AS %s) AS %s", sql_types, q(names(frame))), collapse = ","), "WHERE FALSE")
  DBI::dbExecute(conn, sql_sample_table(conn, name, empty, names(frame)))
  if (!nrow(frame)) return(q(name))
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  batch = max(1L, min(4096L, floor(prediction_state_limit() / ncol(frame))))
  first = 1L
  while (first <= nrow(frame)) {
    last = min(nrow(frame), first + batch - 1L)
    chunk = as.data.frame(lapply(frame[first:last, , drop = FALSE], prediction_frame_text),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (native) {
      values = sql_values_frame(conn, chunk, max_bytes = 8 * 1024^2, format = "json")
      if (is.null(values)) {
        if (batch == 1L) stop("A prediction input row exceeds the SQL transfer limit.", call. = FALSE)
        batch = max(1L, floor(batch / 2L))
        next
      }
      source = paste0("(", values, ") incoming")
      DBI::dbExecute(conn, paste("INSERT INTO", q(name), "SELECT", paste(fields, collapse = ","), "FROM", source))
    } else {
      registered = basename(tempfile("dbreg_prediction_chunk_"))
      tryCatch(
        {
          duckdb::duckdb_register(conn, registered, chunk)
          DBI::dbExecute(conn, paste("INSERT INTO", q(name), "SELECT", paste(fields, collapse = ","), "FROM", q(registered)))
        },
        finally = try(duckdb::duckdb_unregister(conn, registered), silent = TRUE)
      )
    }
    first = last + 1L
  }
  q(name)
}

# Keys are length-prefixed UTF-8 byte strings. Decode only a bounded portable
# dictionary already held in R; no SQL observation data is collected here.
prediction_decode_keys = function(dictionary) {
  width = length(dictionary[["columns"]])
  columns = stats::setNames(rep(list(character(length(dictionary[["key"]]))), width), paste0("k", seq_len(width)))
  for (i in seq_along(dictionary[["key"]])) {
    raw = charToRaw(enc2utf8(dictionary[["key"]][[i]]))
    position = 1L
    for (j in seq_len(width)) {
      remaining = seq.int(position, length(raw))
      colon = remaining[which(raw[remaining] == as.raw(58))[[1L]]]
      size = suppressWarnings(as.numeric(rawToChar(raw[seq.int(position, colon - 1L)])))
      if (!is.finite(size) || size < 0 || size != floor(size) || colon + size > length(raw)) {
        stop("Malformed portable prediction key.", call. = FALSE)
      }
      value = if (size) rawToChar(raw[seq.int(colon + 1L, colon + size)]) else ""
      Encoding(value) = "UTF-8"
      columns[[j]][[i]] = value
      position = colon + size + 1L
    }
    if (position != length(raw) + 1L) stop("Malformed portable prediction key.", call. = FALSE)
  }
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}

prediction_integer_sql_type = function(values, conn) {
  negative = startsWith(values, "-")
  digits = sub("^-", "", values)
  long_limit = ifelse(negative, "9223372036854775808", "9223372036854775807")
  if (all(nchar(digits) < 19L | nchar(digits) == 19L & digits <= long_limit)) return("BIGINT")
  if (max(nchar(digits)) <= 38L) return("DECIMAL(38,0)")
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    stop("Portable integer keys exceed this backend's exact decimal range.", call. = FALSE)
  }
  if (any(negative)) "HUGEINT" else "UHUGEINT"
}
