validate_sample_args = function(row_id, sample_table, conn, table, sql_only, data_only) {
  if (!is.null(row_id) && (!is.character(row_id) || !length(row_id) || anyNA(row_id) || any(!nzchar(row_id)) || anyDuplicated(row_id))) {
    stop("`row_id` must name unique source key columns.", call. = FALSE)
  }
  if (!is.null(sample_table)) {
    if (!is.character(sample_table) || length(sample_table) != 1L || is.na(sample_table) || !nzchar(sample_table)) {
      stop("`sample_table` must name one new temporary table.", call. = FALSE)
    }
    if (is.null(row_id)) stop("`sample_table` requires explicit `row_id` columns.", call. = FALSE)
    if (is.null(conn) && !inherits(table, "tbl_lazy")) stop("`sample_table` requires a caller-owned connection.", call. = FALSE)
    if (isTRUE(sql_only) || isTRUE(data_only)) stop("`sample_table` cannot be combined with sql_only or data_only.", call. = FALSE)
  }
  invisible(NULL)
}

sample_table_exists = function(conn, table) {
  names = if (isTRUE(detect_backend(conn)[["is_databricks"]])) DBI::dbGetQuery(conn, "SHOW TABLES")[["tableName"]] else DBI::dbListTables(conn)
  any(tolower(names) == tolower(table))
}

sql_sample_table = function(conn, table, select_sql, columns) {
  mapping = if (isTRUE(detect_backend(conn)[["is_databricks"]]) && any(grepl("[ ,;{}()\n\r\t=]", columns))) {
    "TBLPROPERTIES ('delta.columnMapping.mode'='name')"
  } else ""
  paste("CREATE TEMP TABLE", DBI::dbQuoteIdentifier(conn, table), mapping, "AS", select_sql)
}

validate_source_row_id = function(conn, from_statement, keys) {
  if (is.null(keys)) return(invisible(NULL))
  backend = detect_backend(conn)[["name"]]
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend))
  if (any(!keys %in% names(schema))) stop("Source row_id columns are unavailable.", call. = FALSE)
  if (folds_identifier_case(conn) && anyDuplicated(tolower(keys))) stop("row_id must name distinct source key columns.", call. = FALSE)
  quoted = as.character(DBI::dbQuoteIdentifier(conn, keys))
  missing = paste0(quoted, " IS NULL")
  numeric = vapply(keys, function(name) is.numeric(schema[[name]]) && !inherits(schema[[name]], "integer64"), logical(1))
  if (backend %in% c("duckdb", "spark") && any(numeric)) missing[numeric] = paste0("(", missing[numeric], " OR isnan(", sql_cast_double(quoted[numeric], conn), "))")
  flag = safe_identifiers("dbreg_invalid_key__", taken = names(schema), fold_case = folds_identifier_case(conn))
  # SQL groups native tuples and tests their counts; no keys or large counts
  # are decoded through R doubles to decide uniqueness.
  sql = paste0("SELECT 1 AS ", flag, " FROM (SELECT * ", from_statement,
    ") dbreg_key_source GROUP BY ", paste(quoted, collapse = ","),
    " HAVING ", sql_count_expr(backend), ">1 OR ", paste(missing, collapse = " OR "))
  if (nrow(DBI::dbGetQuery(conn, sql_limit(sql, 1, backend)))) stop("Source row_id must be unique and nonmissing.", call. = FALSE)
  invisible(NULL)
}

prepare_linear_sample = function(inputs, from_statement, variables) {
  keys = inputs[["row_id"]]
  export = inputs[["sample_table"]]
  if (is.null(keys) && is.null(export)) return(from_statement)
  conn = inputs[["conn"]]
  if (!is.null(export)) {
    backend = detect_backend(conn)
    if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) stop("Sample exports support DuckDB and native Databricks.", call. = FALSE)
    if (isTRUE(inputs[["own_conn"]])) stop("Sample exports require a caller-owned connection.", call. = FALSE)
    if (sample_table_exists(conn, export)) stop("Sample export table already exists.", call. = FALSE)
    repeat {
      snapshot = basename(tempfile("dbreg_sample_"))
      if (tolower(snapshot) != tolower(export) && !sample_table_exists(conn, snapshot)) break
    }
    columns = unique(c(variables, keys))
    sql = paste("SELECT", paste(DBI::dbQuoteIdentifier(conn, columns), collapse = ","), from_statement)
    resources = inputs[["resources"]]
    resources[["tables"]] = c(resources[["tables"]], snapshot)
    DBI::dbExecute(conn, sql_sample_table(conn, snapshot, sql, columns))
    from_statement = paste("FROM", DBI::dbQuoteIdentifier(conn, snapshot))
  }
  validate_source_row_id(conn, from_statement, keys)
  from_statement
}

finalize_linear_sample = function(result, inputs) {
  result[["row_id"]] = inputs[["row_id"]]
  result[["sample_table"]] = inputs[["sample_table"]]
  if (!is.null(inputs[["sample_table"]])) {
    conn = inputs[["conn"]]
    if (sample_table_exists(conn, inputs[["sample_table"]])) stop("Sample export table already exists.", call. = FALSE)
    columns = paste(DBI::dbQuoteIdentifier(conn, inputs[["row_id_safe"]]), "AS", DBI::dbQuoteIdentifier(conn, inputs[["row_id"]]))
    sql = paste("SELECT", paste(columns, collapse = ","), inputs[["from_statement"]])
    # This final CREATE has no replacement clause, and exports never enter
    # the owned-resource tracker. The caller controls their session lifetime.
    DBI::dbExecute(conn, sql_sample_table(conn, inputs[["sample_table"]], sql, inputs[["row_id"]]))
  }
  result
}
