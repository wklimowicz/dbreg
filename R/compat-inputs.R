# Normalize positional R inputs before any native sample filtering.
compat_input_vectors = function(data, table, weights, offset, subset) {
  r_input = is.null(table) && is.data.frame(data)
  if (r_input) data = as.data.frame(data)
  n = if (r_input) nrow(data) else NULL
  add = function(data, value, label) {
    name = safe_identifiers(paste0('dbreg_', label, '__'), taken = names(data), fold_case = TRUE)
    data[[name]] = value
    list(data = data, name = name)
  }
  formula_column = function(value, label) {
    if (!inherits(value, 'formula')) return(value)
    if (length(value) != 2L || !is.symbol(value[[2L]])) stop('`',label,'` formulas must name one column; calculate other expressions in the source.', call. = FALSE)
    as.character(value[[2L]])
  }
  vector_offset = FALSE
  for (label in c('weights','offset')) {
    value = if (label == 'weights') weights else offset
    value = formula_column(value, label)
    if (!is.null(value) && !is.character(value)) {
      if (!r_input) stop('Positional ',label,' vectors require R data; SQL inputs require a column.', call. = FALSE)
      if (!is.numeric(value) || !(length(value) %in% c(1L,n))) stop('`',label,'` must be numeric with length one or the original data size.', call. = FALSE)
      if (length(value) == 1L) value = rep(value, n)
      appended = add(data, value, label)
      data = appended[['data']]
      value = appended[['name']]
      if (label == 'offset') vector_offset = TRUE
    }
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))) {
      stop('`', label, '` must name one column or be a numeric R-data vector.', call. = FALSE)
    }
    if (label == 'weights') weights = value else offset = value
  }
  if (!is.null(subset) && !is.character(subset) && !inherits(subset,'formula')) {
    if (!r_input) stop('Positional subset vectors require R data; SQL inputs require a predicate formula or column.', call. = FALSE)
    if (is.matrix(subset) && ncol(subset) == 1L) subset = as.vector(subset)
    if (is.logical(subset)) {
      if (length(subset) != n) stop('Logical subset vectors must have the original data size.', call. = FALSE)
      keep = subset
    } else if (is.numeric(subset) && length(subset) && all(is.finite(subset)) && all(subset == floor(subset)) && !any(subset == 0)) {
      if (all(subset < 0)) {
        keep = rep(TRUE, n)
        excluded = abs(subset[subset >= -n])
        keep[excluded] = FALSE
      } else if (all(subset > 0) && max(subset) <= n) {
        if (is.unsorted(subset, strictly = TRUE)) stop('Positive subset indices must be unique and in original order; use a logical mask or subset the R data explicitly.', call. = FALSE)
        keep = rep(FALSE, n)
        keep[subset] = TRUE
      } else stop('Subset indices must have one sign and lie within the original data.', call. = FALSE)
    } else stop('Subset must be a logical vector or finite nonzero integer indices.', call. = FALSE)
    appended = add(data, keep, 'subset')
    data = appended[['data']]
    subset = appended[['name']]
  }
  subset_variables(subset)
  list(data = data, weights = weights, offset = offset, subset = subset, vector_offset = vector_offset)
}
