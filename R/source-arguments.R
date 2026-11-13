# Resolve source precedence before touching ignored promises or inspecting data.
normalize_data_source = function(table, data, path) {
  if (!is.null(table)) return(list(table = table, data = NULL, path = NULL))
  if (inherits(data, "tbl_lazy")) return(list(table = data, data = NULL, path = NULL))
  if (!is.null(data)) return(list(table = NULL, data = data, path = NULL))
  list(table = NULL, data = NULL, path = path)
}

formula_column = function(value, label) {
  if (!inherits(value, "formula")) return(value)
  if (length(value) != 2L || !is.symbol(value[[2L]])) {
    stop("`", label, "` formulas must name one column; calculate other expressions in the source.", call. = FALSE)
  }
  as.character(value[[2L]])
}

# Reject misspelled or unsupported statistical requests before touching a
# source or a retained session. Documented legacy names are handled by callers.
check_unused_args = function(dots, allowed = character()) {
  if (!length(dots)) return(invisible(NULL))
  labels = names(dots)
  if (is.null(labels)) labels = rep("", length(dots))
  labels[is.na(labels)] = ""
  unused = !nzchar(labels) | !labels %in% allowed
  if (any(unused)) {
    stop("Unused arguments: ", paste(ifelse(nzchar(labels[unused]),
      paste0("`", labels[unused], "`"), "<unnamed>"), collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}
