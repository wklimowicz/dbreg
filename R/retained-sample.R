retain_prediction = function(retain) length(retain) == 1L && !is.na(retain) && retain %in% c("prediction", "all")

retain_sample = function(retain) length(retain) == 1L && !is.na(retain) && retain %in% c("inference", "all")

prepare_retained_data = function(data, table, retain) {
  r_input = retain_sample(retain) && is.null(table) && is.data.frame(data)
  key = NULL
  if (r_input) {
    data = as.data.frame(data)
    key = safe_identifiers("dbreg_row_key__", taken = names(data), fold_case = TRUE)
    data[[key]] = seq_len(nrow(data))
  }
  list(data = data, key = key, r_input = r_input)
}

materialize_retained_source = function(conn, from_statement, resources, key = NULL, r_input = FALSE) {
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Full fitted-state retention supports DuckDB and native Databricks.", call. = FALSE)
  }
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  schema = DBI::dbGetQuery(conn, paste("SELECT *", from_statement, "LIMIT 0"))
  columns = names(schema)
  generated = is.null(key)
  if (generated) key = safe_identifiers("dbreg_row_key__", taken = columns, fold_case = folds_identifier_case(conn))
  if (!generated && !key %in% columns) stop("The retained source is missing its positional key.", call. = FALSE)
  fields = q(columns)
  if (generated) fields = c(fields, paste(if (isTRUE(backend[["is_databricks"]])) "monotonically_increasing_id()" else
    "ROW_NUMBER() OVER ()", "AS", q(key)))
  name = basename(tempfile("dbreg_original_sample_"))
  resources[["conn"]] = conn
  resources[["tables"]] = c(resources[["tables"]], name)
  DBI::dbExecute(conn, sql_sample_table(conn, name, paste("SELECT", paste(fields, collapse = ","), from_statement), c(columns, if (generated) key)))
  n = as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", q(name)))[["n"]])
  classes = vapply(schema, function(value) class(value)[[1L]], character(1))
  if (generated) classes[[key]] = "integer64"
  list(table = name, key = key, r_input = r_input, columns = setdiff(columns, if (!generated) key), column_classes = classes, nobs_input = n)
}

sql_distinct_tuple_count = function(conn, columns) {
  columns = as.character(DBI::dbQuoteIdentifier(conn, columns))
  if (length(columns) == 1L) return(paste0("COUNT(DISTINCT ", columns, ")"))
  native = isTRUE(detect_backend(conn)[["is_databricks"]])
  paste0("COUNT(DISTINCT ", if (!native) "(" else "",
    paste(columns, collapse = ","), if (!native) ")" else "", ")")
}

# Append scores to the immutable original rows while the fit projector lives.
# Only this final augmented relation is retained, including excluded positions.
retain_score_sample = function(conn, original, resources, score_query, score_key, scores,
  eta, residual, mu = NULL, row_id = NULL, fe_components = list(), weights = NULL) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  private = safe_identifiers(c("dbreg_estimated__", "dbreg_eta__", "dbreg_residual__",
    if (!is.null(mu)) "dbreg_mu__", if (length(scores)) paste0("dbreg_score", seq_along(scores), "__")),
    taken = c(original[["columns"]], original[["key"]]), fold_case = folds_identifier_case(conn))
  estimated_column = private[[1L]]
  eta_column = private[[2L]]
  residual_column = private[[3L]]
  mu_column = if (!is.null(mu)) private[[4L]] else NULL
  score_columns = if (length(scores)) utils::tail(private, length(scores)) else character()
  classes = original[["column_classes"]]
  classes[private] = "numeric"
  classes[[estimated_column]] = "logical"
  select = c("r.*", paste0("s.", q(score_key), " IS NOT NULL AS ", q(estimated_column)),
    paste(eta, "AS", q(eta_column)), paste(residual, "AS", q(residual_column)),
    if (!is.null(mu)) paste(mu, "AS", q(mu_column)), if (length(scores)) paste(scores, "AS", q(score_columns)))
  name = basename(tempfile("dbreg_fitted_sample_"))
  resources[["tables"]] = c(resources[["tables"]], name)
  sql = paste0("SELECT ", paste(select, collapse = ","), " FROM ", q(original[["table"]]), " r LEFT JOIN (",
    score_query, ") s ON r.", q(original[["key"]]), "=s.", q(score_key))
  DBI::dbExecute(conn, sql_sample_table(conn, name, dbglm_sql(conn, sql), c(original[["columns"]], original[["key"]], private)))
  w = if (is.null(weights)) "1.0" else paste0("CAST(", q(weights), " AS DOUBLE)")
  summary = c("COUNT(*) AS n_original", paste0("SUM(CASE WHEN ", q(estimated_column), " THEN 1 ELSE 0 END) AS n_fit"),
    paste0("SUM(CASE WHEN ", q(estimated_column), " THEN POWER(", q(residual_column), ",2) ELSE 0.0 END) AS rss_unweighted"),
    paste0("SUM(CASE WHEN ", q(estimated_column), " THEN ", w, "*POWER(", q(residual_column), ",2) ELSE 0.0 END) AS rss"))
  row = dbglm_get_query(conn, paste("SELECT", paste(summary, collapse = ","), "FROM", q(name)))
  if (as.numeric(row[["n_original"]]) != original[["nobs_input"]]) stop("Retained score rows do not match the original sample.", call. = FALSE)
  fe_sizes = if (length(fe_components)) {
    terms = vapply(seq_along(fe_components), function(i) paste(sql_distinct_tuple_count(conn, fe_components[[i]]), "AS", paste0("g", i)), character(1))
    counts = dbglm_get_query(conn, paste("SELECT", paste(terms, collapse = ","), "FROM", q(name), "WHERE", q(estimated_column)))
    vapply(counts, as.numeric, numeric(1))
  } else numeric()
  list(version = 1L, table = name, key = original[["key"]], r_input = original[["r_input"]],
    row_id = if (length(row_id)) row_id else original[["key"]], estimated = estimated_column,
    eta = eta_column, residual = residual_column, mu = mu_column, scores = score_columns,
    nobs_input = original[["nobs_input"]], nobs = as.numeric(row[["n_fit"]]), rss = as.numeric(row[["rss"]]),
    rss_unweighted = as.numeric(row[["rss_unweighted"]]), fe_components = fe_components, fe_sizes = fe_sizes,
    weights = weights,
    column_classes = classes,
    resources = list(storage = "database", tables = name, version = 1L))
}
