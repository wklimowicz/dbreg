validate_factor_refs = function(refs, variables) {
  if (is.null(refs)) return(NULL)
  refs = as.list(refs)
  if (is.null(names(refs)) || anyNA(names(refs)) || any(!nzchar(names(refs))) ||
      anyDuplicated(names(refs)) || any(!names(refs) %in% variables) ||
      any(lengths(refs) != 1L) || any(vapply(refs, anyNA, logical(1)))) {
    stop("`factor_refs` must name regressor columns with one nonmissing reference each.", call. = FALSE)
  }
  refs
}

factor_reference_index = function(value, keys, type) {
  if (type == "numeric" && is.character(value)) value = suppressWarnings(as.numeric(value))
  if (type == "logical" && !is.logical(value)) {
    value = if (is.numeric(value) && value %in% c(0, 1)) as.logical(value) else
      if (is.character(value) && value %in% c("FALSE", "TRUE")) value == "TRUE" else NA
  }
  position = match(category_key_part(value, type), category_key_part(keys, type))
  if (length(position) != 1L || is.na(position)) {
    stop("Reference level was not found in regressor.", call. = FALSE)
  }
  position
}
