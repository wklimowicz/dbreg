# Parse only explicitly supported one-source transformations. Source names are
# syntax nodes, never strings split on operators or implementation aliases.
parse_linear_term = function(term) {
  expression = str2lang(term)
  if (is.name(expression)) return(list(kind = "numeric", variable = as.character(expression)))
  if (!is.call(expression)) return(NULL)
  if (identical(expression[[1L]], as.name("factor"))) {
    if (length(expression) != 2L || !is.name(expression[[2L]])) {
      stop("Use factor(column) without additional arguments.", call. = FALSE)
    }
    return(list(kind = "factor", variable = as.character(expression[[2L]])))
  }
  if (identical(expression[[1L]], as.name("I"))) {
    power = if (length(expression) == 2L) expression[[2L]] else NULL
    if (!is.call(power) || !identical(power[[1L]], as.name("^")) ||
        length(power) != 3L || !is.name(power[[2L]]) ||
        !is.numeric(power[[3L]]) || length(power[[3L]]) != 1L ||
        is.na(power[[3L]]) || !power[[3L]] %in% c(2, 3)) {
      stop("Only I(column^2) and I(column^3) are supported.", call. = FALSE)
    }
    return(list(kind = "power", variable = as.character(power[[2L]]), exponent = as.integer(power[[3L]])))
  }
  NULL
}

linear_term_key = function(term) {
  paste(deparse(str2lang(term), backtick = TRUE), collapse = "")
}

linear_term_nodes = function(labels) {
  parts = unique(unlist(lapply(labels, split_interaction_term)))
  nodes = lapply(parts, parse_linear_term)
  names(nodes) = vapply(parts, linear_term_key, character(1))
  nodes
}
