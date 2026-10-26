validate_linear_offset = function(offset) {
  if (!is.null(offset) && (!is.character(offset) || length(offset) != 1L || is.na(offset) || !nzchar(offset))) {
    stop("`offset` must name one numeric column.", call. = FALSE)
  }
}

# Derive a distinct outcome column so the original outcome can independently
# occur in regressors, offsets, weights or identity columns.
offset_source_frame = function(conn, from, y, offset) {
  if (is.null(offset)) return(list(from_statement = from))
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("OLS offsets support DuckDB and native Databricks.", call. = FALSE)
  }
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from), 0, backend[["name"]]))
  if (!is.numeric(schema[[offset]]) || !(is.numeric(schema[[y]]) || is.logical(schema[[y]]))) {
    stop("Outcome and offset must be numeric (the outcome may be logical).", call. = FALSE)
  }
  alias = safe_identifiers("dbreg_offset_outcome__", taken = names(schema), fold_case = folds_identifier_case(conn))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  expression = paste0("CAST(", q(y), " AS DOUBLE)-CAST(", q(offset), " AS DOUBLE)")
  from = paste0("FROM (SELECT dbreg_offset_source.*, (", expression, ") AS ", q(alias),
    " FROM (SELECT * ", from, ") dbreg_offset_source) dbreg_offset_frame")
  list(from_statement = from, alias = alias, outcome = y,
    finite = backend_sql(conn, paste0("isfinite(", q(y), ") AND isfinite(", q(offset), ")")))
}

linear_offset_frame = function(conn, from, parsed, offset) {
  if (is.null(offset)) return(list(from_statement = from, parsed = parsed))
  frame = offset_source_frame(conn, from, parsed[["yvar"]], offset)
  formula = stats::as.formula(call("~", as.name(frame[["alias"]]), parsed[["fml"]][[3L]]), env = baseenv())
  frame[["parsed"]] = parse_regression_formula(formula, allow_empty = TRUE, allow_extended = TRUE)
  frame
}

linear_prediction_offset = function(object, data) {
  column = object[["offset"]]
  if (is.null(column)) return(rep(0, nrow(data)))
  if (!column %in% names(data)) stop("newdata is missing required offset column: ", column, call. = FALSE)
  value = data[[column]]
  if (!is.numeric(value)) stop("The prediction offset must be numeric.", call. = FALSE)
  value = as.numeric(value)
  value[!is.finite(value)] = NA_real_
  value
}
