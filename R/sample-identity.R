# Source keys are typed tuples. Never infer their alignment from fetch order.
validate_row_id_arg = function(row_id) {
  if (!is.null(row_id) && (!is.character(row_id) || !length(row_id) ||
    anyNA(row_id) || any(!nzchar(row_id)) || anyDuplicated(row_id))) {
    stop("`row_id` must name unique source key columns.", call. = FALSE)
  }
  invisible(NULL)
}

validate_source_row_id = function(conn, from_statement, row_id) {
  if (is.null(row_id)) return(invisible(NULL))
  quoted = as.character(DBI::dbQuoteIdentifier(conn, row_id))
  counter = safe_identifiers("dbreg_key_count__", taken = row_id,
    fold_case = folds_identifier_case(conn))
  count = sql_count_expr(detect_backend(conn)[["name"]])
  missing = paste0(quoted, " IS NULL")
  backend = detect_backend(conn)
  if (backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) {
    text = if (isTRUE(backend[["is_databricks"]])) "STRING" else "VARCHAR"
    missing = c(missing, paste0("(LOWER(typeof(", quoted, ")) IN ('double','float','real') AND LOWER(CAST(",
      quoted, " AS ", text, "))='nan')"))
  }
  # Compare counts in SQL so even counts outside double's exact integer range
  # cannot turn a duplicate into a false successful validation in R.
  info = DBI::dbGetQuery(conn, paste0("SELECT CASE WHEN COALESCE(SUM(", counter,
    "),0)=", count, " THEN 0 ELSE 1 END AS duplicate_keys,COALESCE(SUM(CASE WHEN ",
    paste(missing, collapse = " OR "), " THEN ", counter,
    " ELSE 0 END),0) AS n_null FROM (SELECT ", paste(quoted, collapse = ","),
    ",", sql_count(conn, counter), " ", from_statement,
    " GROUP BY ", paste(quoted, collapse = ","), ") dbreg_keys"))
  if (info[["n_null"]] > 0 || info[["duplicate_keys"]] > 0) {
    stop("Source row_id must be unique and nonmissing.", call. = FALSE)
  }
  invisible(NULL)
}

# Both supported export backends accept this CREATE form. Deliberately avoid
# OR REPLACE: an export must never overwrite a caller-owned relation.
sample_table_exists = function(conn, table) {
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    # ODBC's SQLTables metadata can omit session tables. Do not use SHOW LIKE:
    # table names may contain literal pattern metacharacters such as * or |.
    tables = DBI::dbGetQuery(conn, "SHOW TABLES")[["tableName"]]
    return(any(tolower(tables) == tolower(table)))
  }
  DBI::dbExistsTable(conn, table)
}

sql_sample_table = function(conn, table, select_sql, columns, replace = FALSE) {
  # Delta's default Parquet naming constraints also apply to session tables.
  # Enable mapping only on our new tables when preserving such source names.
  # https://docs.databricks.com/aws/en/tables/features/column-mapping
  mapping = if (isTRUE(detect_backend(conn)[["is_databricks"]]) &&
    any(grepl("[ ,;{}()\n\r\t=]", columns))) {
    "TBLPROPERTIES ('delta.columnMapping.mode'='name')"
  } else ""
  paste("CREATE", if (replace) "OR REPLACE" else "", "TEMP TABLE",
    DBI::dbQuoteIdentifier(conn, table), mapping, "AS", select_sql)
}

export_retained_sample = function(conn, table, keys, from_statement, public_names = keys) {
  columns = paste(DBI::dbQuoteIdentifier(conn, keys), "AS",
    DBI::dbQuoteIdentifier(conn, public_names))
  sql = paste("SELECT", paste(columns, collapse = ","), from_statement)
  DBI::dbExecute(conn, sql_sample_table(conn, table, sql, public_names))
  invisible(NULL)
}

# Explicit sample exports use an owned snapshot for the whole fit. Default
# compact fits keep their existing query and resource behavior.
prepare_linear_sample = function(inputs, from_statement, variables) {
  keys = inputs[["row_id"]]
  table = inputs[["sample_table"]]
  full = retain_sample(inputs[["retain"]])
  prediction = retain_prediction(inputs[["retain"]]) && length(inputs[["fe_labels"]]) > 0L
  if (is.null(keys) && is.null(table) && !prediction && !full) return(from_statement)
  conn = inputs[["conn"]]
  if (full && !is.null(table) && sample_table_exists(conn, table)) stop("Sample export table already exists.", call. = FALSE)
  if (full) {
    original = materialize_retained_source(conn, from_statement, inputs[["sample_resources"]],
      inputs[["retained_input"]][["key"]], inputs[["retained_input"]][["r_input"]])
    inputs[["retained_source"]] = original
    from_statement = paste("FROM", DBI::dbQuoteIdentifier(conn, original[["table"]]))
  }
  if ((!is.null(table) || prediction) && !full) {
    backend = detect_backend(conn)
    if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
      stop("Sample snapshots currently support DuckDB and native Databricks.", call. = FALSE)
    }
    if (!is.null(table) && sample_table_exists(conn, table)) stop("Sample export table already exists.", call. = FALSE)
    resources = inputs[["sample_resources"]]
    resources[["conn"]] = conn
    columns = unique(c(variables, keys))
    snapshot = basename(tempfile("dbreg_sample_"))
    # Track before CREATE so partially failed setup is cleaned by on.exit too.
    resources[["tables"]] = c(resources[["tables"]], snapshot)
    sql = paste("SELECT", paste(DBI::dbQuoteIdentifier(conn, columns), collapse = ","), from_statement)
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
    export_retained_sample(inputs[["conn"]], inputs[["sample_table"]],
      inputs[["row_id_safe"]], inputs[["from_statement"]], inputs[["row_id"]])
  }
  result
}
