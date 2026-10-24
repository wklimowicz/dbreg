# Materialized, explicitly keyed output. No input ordering is assumed and no
# prediction panel is implicitly collected into R.
predict_sql = function(object, newdata, conn, table, path, row_id, prediction_table, type = "link") {
  validate_row_id_arg(row_id)
  if (is.null(row_id)) stop("Database prediction requires explicit `row_id` columns.", call. = FALSE)
  if (inherits(object, "dbivreg")) stop("Database prediction output is not yet supported for IV models.", call. = FALSE)
  if (identical(object[["strategy"]], "mundlak")) stop("Database prediction output is not yet supported for Mundlak models.", call. = FALSE)
  if (sum(!vapply(list(newdata, table, path), is.null, logical(1))) != 1L) {
    stop("Provide exactly one of `newdata`, `table`, or `path` for database prediction.", call. = FALSE)
  }
  if (!is.null(prediction_table) && (!is.character(prediction_table) || length(prediction_table) != 1L ||
    is.na(prediction_table) || !nzchar(prediction_table))) stop("`prediction_table` must be a single nonempty table name.", call. = FALSE)
  caller_connection = !is.null(conn)
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) fitted_resource_connection(inherited)
  state = object[["fitted_state"]]
  database = identical(state[["prediction"]][["storage"]], "database")
  if (database) {
    owner = fitted_resource_connection(object[["resource_state"]])
    if (!is.null(conn) && !fitted_resource_same_connection(conn, owner)) {
      stop("Database-backed prediction requires its fitted owning session.", call. = FALSE)
    }
    conn = owner
  }
  if (!is.null(prediction_table) && (!caller_connection ||
    (!is.null(inherited) && isTRUE(inherited[["session"]][["own_conn"]])) ||
    (database && isTRUE(object[["resource_state"]][["session"]][["own_conn"]])))) {
    stop("An explicit prediction export requires a caller-owned connection.", call. = FALSE)
  }
  required = prediction_required_columns(object)
  labels = NULL
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) stop("`newdata` must be a data frame.", call. = FALSE)
    incoming = prediction_dataframe_source(object, newdata, unique(c(required, row_id)))
    newdata = incoming[["data"]]
    labels = incoming[["labels"]]
  }
  resources = dbglm_resources()
  on.exit(dbglm_release_resources(resources), add = TRUE)
  if (database) resources[["retained_session"]] = object[["resource_state"]][["session"]]
  # Upload incoming R values as typed text on both engines. This also avoids
  # bigint="numeric" registrations reinterpreting integer64 storage bits.
  if (!is.null(newdata)) {
    if (is.null(conn)) {
      conn = DBI::dbConnect(duckdb::duckdb(bigint = "integer64"), shutdown = TRUE)
      resources[["own_conn"]] = TRUE
    }
    resources[["conn"]] = conn
    inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
    if (!is.null(inherited)) fitted_resource_connection(inherited)
    from_statement = paste("FROM", prediction_upload_frame(conn, newdata, resources))
  } else {
    if (is.null(conn) && !is.null(table) && !inherits(table, "tbl_lazy")) stop("A table prediction source requires `conn`.", call. = FALSE)
    setup = setup_db_connection(conn, table, NULL, path, caller = "prediction", resources = resources)
    conn = setup[["conn"]]
    from_statement = setup[["from_statement"]]
  }
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Database prediction output currently supports DuckDB and native Databricks.", call. = FALSE)
  }
  if (!is.null(prediction_table) && sample_table_exists(conn, prediction_table)) stop("Prediction export table already exists.", call. = FALSE)
  context = sql_prediction_context(object, conn, labels)
  source = sql_prediction_source(object, context, conn, from_statement, row_id, required, labels)
  dictionaries = prediction_portable_tables(context, conn, resources)
  prediction = safe_identifiers("prediction", taken = row_id, fold_case = folds_identifier_case(conn))
  query = sql_prediction_equation(object, context, source, dictionaries, conn, row_id, prediction)
  create = function(query) {
    name = basename(tempfile("dbreg_prediction_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, sql_sample_table(conn, name, dbglm_sql(conn, query), c(row_id, prediction)))
    name
  }
  output = create(query)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  validate_source_row_id(conn, paste("FROM", q(output)), row_id)
  if (inherits(object, "dbglm") && type == "response") {
    keys = paste0("dbreg_prediction_key", seq_along(row_id), "__")
    mean_source = paste0("(SELECT ", paste(q(row_id), "AS", keys, collapse = ","), ",", q(prediction),
      " AS __prediction_eta,CAST(0.0 AS DOUBLE) AS __dbglm_y,CAST(1.0 AS DOUBLE) AS __dbglm_a,",
      "CAST(0.0 AS DOUBLE) AS __dbglm_o FROM ", q(output), ")")
    working = dbglm_working_sql(mean_source, character(), numeric(), object[["family"]][["link"]], eta = "__prediction_eta")
    output = create(paste0(working, "SELECT ", paste(keys, "AS", q(row_id), collapse = ","),
      ",CASE WHEN __prediction_eta IS NULL THEN CAST(NULL AS DOUBLE) ELSE mu END AS ", q(prediction), " FROM means"))
  }
  if (is.null(prediction_table)) {
    lease = retain_fit_resources(resources, output)
  } else {
    DBI::dbExecute(conn, sql_sample_table(conn, prediction_table, paste("SELECT * FROM", q(output)), c(row_id, prediction)))
    output = prediction_table
    lease = retain_fit_resources(resources, character())
  }
  structure(list(conn = fitted_resource_public_connection(lease), table = output, row_id = row_id,
    prediction_column = prediction, resource_state = lease), class = "dbreg_prediction")
}

#' Database prediction results
#'
#' Keyed prediction output stays in its database session. `head()` collects a
#' bounded preview; `as.data.frame()` explicitly collects the complete result.
#' `close()` releases owned output tables. Caller-named exports survive release.
#' Output leases share their owning session, so closing a fitted model does not
#' invalidate its prediction results. A model-owned session closes after its
#' final lease is released. Live results cannot be restored by serialization.
#'
#' @param x,con A `dbreg_prediction` result.
#' @param n Number of rows to preview.
#' @param ... Additional arguments are not supported.
#' @return `as.data.frame()` and `head()` return a data frame containing source
#'   keys and the prediction column. `close()` returns `NULL` invisibly.
#' @name dbreg-prediction
NULL

collect_sql_prediction = function(x, n = NULL) {
  conn = fitted_resource_connection(x[["resource_state"]])
  q = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  from = paste("FROM", q(x[["table"]]), if (!is.null(n)) paste("LIMIT", format(n, scientific = FALSE)) else "")
  output = DBI::dbGetQuery(conn, paste("SELECT *", from))
  if (!nrow(output)) return(output)
  keys = x[["row_id"]]
  types = DBI::dbGetQuery(conn, paste("SELECT", paste0("typeof(", q(keys), ") AS t", seq_along(keys), collapse = ","), "FROM", q(x[["table"]]), "LIMIT 1"))
  exact_text = vapply(seq_along(keys), function(i) {
    type = tolower(types[[i]][[1L]])
    value = output[[keys[[i]]]]
    startsWith(type, "decimal") || (prediction_sql_kind(type) == "integer" &&
      is.numeric(value) && !inherits(value, "integer64") && any(abs(value) >= 2^53, na.rm = TRUE))
  }, logical(1))
  if (any(exact_text)) {
    fields = names(output)
    text = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "STRING" else "VARCHAR"
    select = q(fields)
    index = match(keys[exact_text], fields)
    select[index] = paste0("CAST(", q(fields[index]), " AS ", text, ") AS ", q(fields[index]))
    output = DBI::dbGetQuery(conn, paste("SELECT", paste(select, collapse = ","), from))
  }
  output
}

#' @rdname dbreg-prediction
#' @export
as.data.frame.dbreg_prediction = function(x, ...) {
  check_unused_args(list(...))
  collect_sql_prediction(x)
}

#' @rdname dbreg-prediction
#' @importFrom utils head
#' @export
head.dbreg_prediction = function(x, n = 6L, ...) {
  check_unused_args(list(...))
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) stop("`n` must be a nonnegative whole number.", call. = FALSE)
  collect_sql_prediction(x, n)
}

#' @rdname dbreg-prediction
#' @export
print.dbreg_prediction = function(x, ...) {
  check_unused_args(list(...))
  cat("Database predictions in ", x[["table"]], "; keys: ", paste(x[["row_id"]], collapse = ", "),
    "; prediction column: ", x[["prediction_column"]], ".\n", sep = "")
  invisible(x)
}

#' @rdname dbreg-prediction
#' @export
close.dbreg_prediction = function(con, ...) {
  check_unused_args(list(...))
  fitted_resource_release(con[["resource_state"]])
}
