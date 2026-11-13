# Cohort projections use exact group identities and normalized nonnegative
# weights. They do not recover or store training-sample fixed effects.
prediction_group_index = function(value) {
  if (is.factor(value)) value = as.character(value)
  if (inherits(value, "POSIXt")) value = as.numeric(value)
  type = if (inherits(value, "integer64") || is.integer(value)) "integer" else
    if (inherits(value, "Date")) "date" else if (is.logical(value)) "logical" else
      if (is.numeric(value)) "numeric" else "character"
  keys = category_key_part(value, type)
  if (anyNA(keys)) stop("Prediction grouping values must be finite and nonmissing.", call. = FALSE)
  match(keys, unique(keys))
}

prediction_weights = function(data, name) {
  if (is.null(name)) return(rep(1, nrow(data)))
  value = data[[name]]
  if (!(is.numeric(value) || is.logical(value)) || !is.null(dim(value))) {
    stop("Prediction weights must be numeric.", call. = FALSE)
  }
  value = as.numeric(value)
  if (any(!is.finite(value) | value < 0)) stop("Prediction weights must be finite and non-negative.", call. = FALSE)
  if (!length(value)) return(numeric())
  scale = max(value)
  if (scale <= 0) return(rep(0, length(value)))
  value / scale
}

weighted_group_mean = function(x, group, weights = NULL) {
  if (is.null(weights)) weights = rep(1, length(x))
  weights = as.numeric(weights)
  out = rep(NA_real_, length(x))
  for (rows in split(seq_along(x), prediction_group_index(group))) {
    w = weights[rows]
    if (!length(w) || max(w) <= 0) next
    w = w / max(w)
    w = w / sum(w)
    values = x[rows]
    anchor = min(values)/2 + max(values)/2
    out[rows] = anchor + sum(w * (values - anchor))
  }
  out
}

prepare_prediction_cohort = function(object, data) {
  required = unique(c(model_vars(stats::formula(object[["fml"]], lhs = 0, rhs = 1)),
    if (is.null(object[["fe_components"]])) object[["fe"]] else unlist(object[["fe_components"]], use.names = FALSE), object[["weights"]], object[["offset"]],
    if (object[["strategy"]] == "demean" && object[["yvar"]] %in% names(data)) object[["yvar"]]))
  missing = setdiff(required, names(data))
  if (length(missing)) stop("newdata is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  valid = rep(TRUE, nrow(data))
  for (name in required) {
    value = data[[name]]
    valid = valid & !is.na(value)
    if (is.numeric(value) && !inherits(value, "integer64")) valid = valid & is.finite(value)
  }
  rows = which(valid)
  selected = data[rows, , drop = FALSE]
  if (length(rows) && !is.null(object[["linear_design"]])) {
    keep = !linear_category_invalid(object[["linear_design"]], selected)
    rows = rows[keep]
    selected = selected[keep, , drop = FALSE]
  }
  if (length(rows)) {
    weights = prediction_weights(selected, object[["weights"]])
    supported = rep(TRUE, nrow(selected))
    for (value in prediction_fe_groups(object, selected)) {
      group = prediction_group_index(value)
      totals = rowsum(weights, group, reorder = FALSE)
      supported = supported & totals[match(group, as.integer(rownames(totals))), 1L] > 0
    }
    rows = rows[supported]
  }
  valid[] = FALSE
  valid[rows] = TRUE
  list(data = data[rows, , drop = FALSE], valid = valid)
}

# Alternating weighted projections handle unbalanced and multi-way cohorts.
# Certify every group's residual mean relative to each projected column's RMS.
demean_prediction_data = function(value, groups, weights = NULL, tol = 1e-12, max_iter = 10000L) {
  value = as.matrix(value)
  if (!length(groups) || !nrow(value) || !ncol(value)) return(value)
  if (any(!is.finite(value))) stop("Prediction projection requires finite design values.", call. = FALSE)
  if (is.null(weights)) weights = rep(1, nrow(value))
  weights = as.numeric(weights)
  weights = weights / max(weights)
  centers = vapply(seq_len(ncol(value)), function(j) weighted_group_mean(value[,j], rep(1, nrow(value)), weights)[1L], numeric(1))
  value = sweep(value, 2L, centers, "-")
  scales = apply(abs(value), 2L, max)
  scales[scales == 0] = 1
  value = sweep(value, 2L, scales, "/")
  indices = lapply(groups, prediction_group_index)
  totals = lapply(indices, function(group) rowsum(weights, group, reorder = FALSE))
  means = function(matrix, k) {
    summed = rowsum(matrix * weights, indices[[k]], reorder = FALSE)
    averaged = summed / as.numeric(totals[[k]])
    averaged[match(indices[[k]], as.integer(rownames(summed))), , drop = FALSE]
  }
  for (iteration in seq_len(max_iter)) {
    for (k in seq_along(indices)) value = value - means(value, k)
    rms = sqrt(colSums(value^2 * weights) / sum(weights))
    rms[rms == 0] = 1
    certificate = 0
    for (k in seq_along(indices)) certificate = max(certificate, abs(sweep(means(value, k), 2L, rms, "/")))
    if (is.finite(certificate) && certificate <= tol) return(sweep(value, 2L, scales, "*"))
  }
  stop("Prediction demeaning did not converge.", call. = FALSE)
}
