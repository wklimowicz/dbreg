validate_diagnostics = function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}

# Optional bounded aggregate over the original source. Counts overlap across
# roles and conditions; only the fitted count is the final estimation sample.
sample_diagnostics = function(inputs, from_statement, variables, source) {
  if (!isTRUE(inputs[["diagnostics"]]) || isTRUE(inputs[["sql_only"]]) || isTRUE(inputs[["data_only"]])) return(NULL)
  conn = inputs[["conn"]]
  backend = detect_backend(conn)[["name"]]
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend))
  roles = rep(names(variables), lengths(variables))
  columns = unlist(variables, use.names = FALSE)
  if (any(!columns %in% names(schema))) stop("Audit source columns are unavailable.", call. = FALSE)
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  numeric = vapply(columns, function(name) is.numeric(schema[[name]]) && !inherits(schema[[name]], "integer64"), logical(1))
  missing = paste0(q(columns), " IS NULL")
  if (backend %in% c("duckdb", "spark") && any(numeric)) {
    missing[numeric] = paste0("(", missing[numeric], " OR isnan(", sql_cast_double(q(columns[numeric]), conn), "))")
  }
  # Allocate output aliases separately from every source role. CASE expressions
  # cannot accidentally bind to an earlier aggregate alias with the same name.
  labels = c("n_input", "n_selected", "n_nonpositive_weight", "n_invalid_offset", paste0("n_missing_", seq_along(columns)))
  aliases = setNames(safe_identifiers(paste0("dbreg_audit_", labels), taken = names(schema), fold_case = folds_identifier_case(conn)), labels)
  count_when = function(condition, name) sql_count(conn, aliases[[name]], paste0("CASE WHEN (", condition, ") THEN 1 END"))
  subset = sql_subset_predicate(conn, inputs[["subset"]], schema, source[["to_safe"]])
  if (is.null(subset)) subset = "1=1"
  terms = c(sql_count(conn, aliases[["n_input"]]), count_when(subset, "n_selected"))
  weights = variables[["weights"]]
  terms = c(terms, count_when(if (length(weights)) paste0(q(weights), " <= 0") else "1=0", "n_nonpositive_weight"))
  offset = variables[["offset"]]
  invalid_offset = "1=0"
  if (length(offset)) {
    outcome = variables[["outcome"]]
    values = c(sql_cast_double(q(outcome), conn), sql_cast_double(q(offset), conn),
      paste0("(", sql_cast_double(q(outcome), conn), "-", sql_cast_double(q(offset), conn), ")"))
    finite = paste0("ABS(", values, ") < CAST('Infinity' AS DOUBLE)")
    invalid_offset = paste0("NOT COALESCE((", paste(finite, collapse = " AND "), "),FALSE)")
  }
  terms = c(terms, count_when(invalid_offset, "n_invalid_offset"),
    vapply(seq_along(columns), function(i) count_when(missing[[i]], paste0("n_missing_", i)), character(1)))
  sql = paste0("SELECT ", paste(terms, collapse = ","), " FROM (SELECT * ", from_statement, ") dbreg_audit_source")
  counts = DBI::dbGetQuery(conn, sql)
  read_count = function(name) {
    value = as.numeric(counts[[aliases[[name]]]])
    if (length(value) != 1L || !is.finite(value) || value < 0) stop("Invalid sample audit count.", call. = FALSE)
    value
  }
  list(n_input = read_count("n_input"), n_selected = read_count("n_selected"),
    n_nonpositive_weight = read_count("n_nonpositive_weight"), n_invalid_offset = read_count("n_invalid_offset"),
    detail = data.frame(role = roles, variable = source[["to_display"]](columns),
      n_missing = vapply(paste0("n_missing_", seq_along(columns)), read_count, numeric(1)), row.names = NULL),
    query_string = sql, count_aliases = aliases)
}

finalize_sample_diagnostics = function(result, inputs) {
  audit = inputs[["sample_diagnostics"]]
  result[["nobs_input"]] = if (is.null(audit)) NA_real_ else audit[["n_input"]]
  if (!is.null(audit)) {
    audit[["n_fitted"]] = as.numeric(result[["nobs_orig"]])
    audit[["n_dropped"]] = audit[["n_input"]] - audit[["n_fitted"]]
    result[["sample_diagnostics"]] = audit
  }
  result
}
