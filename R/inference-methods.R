inference_table = function(estimate, covariance, df, normal = FALSE) {
  se = sqrt(Matrix::diag(covariance))
  statistic = estimate / se
  p = if (normal) 2 * stats::pnorm(-abs(statistic)) else 2 * stats::pt(-abs(statistic), df)
  cbind(estimate = estimate, std.error = se, statistic = statistic, p.value = p)
}

normal_inference = function(object) {
  inherits(object, "dbglm") && object[["family"]][["family"]] %in% c("binomial", "poisson")
}

inference_labels = function(object) {
  if (normal_inference(object)) c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  else c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
}

coefficient_view = function(object, fe, collin = FALSE) {
  layout = object[["linear_design"]]
  columns = layout[["columns"]]
  index = if (collin) seq_len(nrow(columns)) else match(layout[["fitted_ids"]], columns[["id"]])
  if (!isTRUE(fe) && length(object[["fe"]])) {
    index = index[columns[["role"]][index] == "slope" & columns[["alias"]][index] != "(Intercept)"]
  }
  list(labels = columns[["label"]][index], fitted = match(columns[["id"]][index], layout[["fitted_ids"]]))
}

coefficient_limits = function(table, level, df, normal) {
  validate_confidence_level(level)
  tail = (1 - level) / 2
  critical = if (normal) stats::qnorm(1 - tail) else stats::qt(1 - tail, df)
  half = critical * table[, "std.error"]
  cbind(table[, "estimate"] - half, table[, "estimate"] + half)
}

visible_coefficient_table = function(object, fe = FALSE) {
  table = object[["coeftable"]]
  if (isTRUE(fe) || !length(object[["fe"]])) return(table)
  positions = if (!is.null(object[["linear_design"]])) visible_linear_positions(object, fe)
    else names(coef(object, fe = fe))
  table[positions, , drop = FALSE]
}

#' Summarize a database regression
#' @param object A fitted `dbreg` or `dbglm` model.
#' @param vcov,cluster Covariance type or cluster terms as in \link{vcov.dbreg}.
#' @param ssc Optional small-sample correction policy. A list may specify
#'   `K.adj`, `K.fixef`, `G.adj`, `G.df` and `t.df`.
#' @param vcov_fix Repair an indefinite two-way covariance matrix.
#' @param ... Unused.
#' @return A fitted-model summary with the selected covariance and coefficient
#'   table. A new covariance request uses the retained estimation sample.
#' @export
summary.dbreg = function(object, vcov = NULL, cluster = NULL, ssc = NULL, vcov_fix = TRUE, ...) {
  check_unused_args(list(...))
  requested = !is.null(vcov) || !is.null(cluster) || !is.null(ssc) || !identical(vcov_fix, TRUE)
  if (requested) {
    V = vcov.dbreg(object, fe = TRUE, vcov = vcov, cluster = cluster, ssc = ssc, vcov_fix = vcov_fix)
    estimate = coef(object, fe = TRUE)
    object[["coeftable"]] = inference_table(estimate, V, attr(V, "df.t"), normal_inference(object))
    object[["vcov"]] = V
    prediction = attr(V, "prediction_covariance", exact = TRUE)
    if (!is.null(object[["prediction_coordinates"]]) && !is.null(prediction)) {
      object[["prediction_coordinates"]][["covariance"]] = prediction
    }
    if (!is.null(vcov) || !is.null(cluster)) {
      if (inherits(vcov, "formula")) environment(vcov) = baseenv()
      if (inherits(cluster, "formula")) environment(cluster) = baseenv()
      object[["postfit_request"]] = list(vcov = vcov, cluster = cluster)
    }
  }
  class(object) = unique(c("summary.dbreg", class(object)))
  object
}

#' Coefficient summaries through fixest generics
#' @param object A fitted database regression.
#' @param keep,drop,order Optional regular expression selectors for coefficient
#'   names. A leading `!` reverses a match.
#' @param ... Covariance arguments passed to \link{summary.dbreg}.
#' @return A coefficient table or a named numeric vector.
#' @name coefficient-summaries
NULL

#' @rdname coefficient-summaries
#' @exportS3Method fixest::coeftable
coeftable.dbreg = function(object, keep = NULL, drop = NULL, order = NULL, ...) {
  object = summary.dbreg(object, ...)
  table = visible_coefficient_table(object)
  labels = rownames(table)
  positions = seq_len(nrow(table))
  match = function(pattern) {
    pattern = sub("^!%", "%!", pattern)
    pattern = sub("^%", "", pattern)
    invert = startsWith(pattern, "!")
    selected = grepl(if (invert) substring(pattern, 2L) else pattern, labels, perl = TRUE)
    if (invert) !selected else selected
  }
  for (value in list(keep, drop, order)) {
    if (!is.null(value) && (!is.character(value) || anyNA(value))) {
      stop("Coefficient selectors must be nonmissing regular expressions.", call. = FALSE)
    }
  }
  if (length(keep)) positions = positions[Reduce(`|`, lapply(keep, match))]
  if (length(drop)) positions = positions[!Reduce(`|`, lapply(drop, match))[positions]]
  for (pattern in rev(order)) {
    first = match(pattern)[positions]
    positions = c(positions[first], positions[!first])
  }
  table = as.matrix(table[positions, , drop = FALSE])
  colnames(table) = inference_labels(object)
  table
}

#' @rdname coefficient-summaries
#' @exportS3Method fixest::se
se.dbreg = function(object, ...) coefficient_statistic(object, 2L, ...)

#' @rdname coefficient-summaries
#' @exportS3Method fixest::tstat
tstat.dbreg = function(object, ...) coefficient_statistic(object, 3L, ...)

#' @rdname coefficient-summaries
#' @exportS3Method fixest::pvalue
pvalue.dbreg = function(object, ...) coefficient_statistic(object, 4L, ...)

coefficient_statistic = function(object, column, ...) {
  table = coeftable.dbreg(object, ...)
  setNames(as.numeric(table[, column]), rownames(table))
}

validate_confidence_level = function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must lie strictly between zero and one.", call. = FALSE)
  }
  invisible(NULL)
}
