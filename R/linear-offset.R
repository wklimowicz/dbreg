linear_offset_column = function(value) {
  if (is.null(value)) return(NULL)
  value = formula_column(value, "offset")
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("`offset` must name one numeric column.", call. = FALSE)
  }
  value
}

# Preserve raw outcome/offset roles independently of the adjusted target.
prepare_linear_offset = function(inputs, from_statement, yvar, source) {
  if (is.null(inputs[["offset"]])) return(list(from_statement = from_statement, yvar = yvar))
  conn = inputs[["conn"]]
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("OLS/IV offsets support DuckDB and native Databricks.", call. = FALSE)
  }
  offset = source[["to_safe"]](inputs[["offset"]])
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend[["name"]]))
  if (!offset %in% names(schema)) stop("Offset column not found: ", inputs[["offset"]], call. = FALSE)
  if (!is.numeric(schema[[offset]]) || !(is.numeric(schema[[yvar]]) || is.logical(schema[[yvar]]))) {
    stop("Outcome and offset must be numeric (the outcome may be logical).", call. = FALSE)
  }
  alias = safe_identifiers("dbreg_offset_outcome__", taken = names(schema), fold_case = folds_identifier_case(conn))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  expression = paste0(sql_cast_double(q(yvar), conn), "-", sql_cast_double(q(offset), conn))
  from_statement = paste0("FROM (SELECT dbreg_offset_source.*, (", expression, ") AS ", q(alias),
    " FROM (SELECT * ", from_statement, ") dbreg_offset_source) dbreg_offset_frame")
  finite = paste0("ABS(", sql_cast_double(q(c(yvar, offset, alias)), conn), ") < CAST('Infinity' AS DOUBLE)")
  from_statement = filter_source_predicate(from_statement, paste(finite, collapse = " AND "))
  inputs[["offset_outcome"]] = source[["to_display"]](yvar)
  inputs[["offset_safe"]] = offset
  list(from_statement = from_statement, yvar = alias)
}

linear_prediction_offset = function(object, data) {
  name = object[["offset"]]
  if (is.null(name)) return(rep(0, nrow(data)))
  if (!name %in% names(data)) stop("newdata is missing required offset column: ", name, call. = FALSE)
  value = data[[name]]
  if (!is.numeric(value)) stop("The prediction offset must be numeric.", call. = FALSE)
  value = as.numeric(value)
  value[!is.finite(value)] = NA_real_
  value
}
