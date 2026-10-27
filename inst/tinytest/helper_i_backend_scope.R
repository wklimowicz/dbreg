# Existing typed SQL backends need no driver adapter. This fixture evaluates
# compiler SQL in DuckDB while exposing declared server type names at typeof().
# It proves dispatch/schema/SQL/coding contracts, not remote engine execution.
dbreg_test_existing_i_backends = function() {
  connection_class = "dbreg_i_existing_backend_fixture"
  if (!methods::isClass(connection_class)) methods::setClass(connection_class,
    contains = "DBIConnection", slots = c(native = "ANY", state = "environment"))
  methods::setMethod("dbGetInfo", connection_class, function(dbObj, ...) list(dbms.name = dbObj@state[["server"]]))
  methods::setMethod("dbQuoteIdentifier", c(connection_class, "character"), function(conn, x, ...) {
    if (identical(conn@state[["server"]], "BigQuery")) DBI::SQL(paste0("`", x, "`")) else DBI::dbQuoteIdentifier(conn@native, x, ...)
  })
  methods::setMethod("dbQuoteLiteral", c(connection_class, "character"), function(conn, x, ...) DBI::dbQuoteLiteral(conn@native, x, ...))
  methods::setMethod("dbQuoteLiteral", c(connection_class, "Date"), function(conn, x, ...) DBI::dbQuoteLiteral(conn@native, x, ...))
  methods::setMethod("dbGetQuery", c(connection_class, "character"), function(conn, statement, ...) {
    conn@state[["queries"]] = c(conn@state[["queries"]], statement)
    if (startsWith(statement, "SELECT typeof(")) {
      expected = paste0("SELECT typeof(", DBI::dbQuoteIdentifier(conn, names(conn@state[["types"]])), ") AS t1 FROM ")
      position = which(startsWith(statement, expected))
      if (length(position) != 1L) stop("Unexpected metadata query: ", statement)
      return(data.frame(t1 = unname(conn@state[["types"]][position])))
    }
    # Only the fixture translates recorded BigQuery identifier spelling.
    statement = gsub("`([^`]+)`", '\"\\1\"', statement)
    DBI::dbGetQuery(conn@native, statement, ...)
  })
  methods::setMethod("dbSendQuery", c(connection_class, "character"), function(conn, statement, ...) {
    stop("Existing typed backend must retain its dbGetQuery schema route.")
  })
  methods::setMethod("dbExecute", c(connection_class, "character"), function(conn, statement, ...) {
    stop("Compiling i() must not acquire a physical resource.")
  })
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  old = options(dbreg.i_max_columns = 1000L)
  on.exit(options(old), add = TRUE)
  data = data.frame(id = seq_len(12), text = rep(c("Beta", "alpha", "gamma"), 4),
    raw_key = rep(c("-9223372036854775808", "9007199254740993", "9223372036854775807"), 4),
    adjacent = rep(1 + c(0, 1, 2) * .Machine[["double.eps"]], 4),
    flag = rep(c(FALSE, TRUE), 6), day = as.Date("2020-01-01") + rep(0:2, 4))
  DBI::dbWriteTable(native, "existing_backend_source", data)
  DBI::dbExecute(native, "CREATE TABLE existing_backend_input AS SELECT *, CAST(raw_key AS BIGINT) AS key FROM existing_backend_source")
  state = new.env(parent = emptyenv())
  state[["queries"]] = character()
  state[["types"]] = c(text = "varchar(8)", key = "bigint", adjacent = "double", flag = "boolean", day = "date")
  conn = methods::new(connection_class, native = native, state = state)
  # Simulate the declared BigQuery connection-class identity through the actual
  # metadata owner, without defining or claiming a real external driver class.
  scope = new.env(parent = asNamespace("dbreg"))
  scope[["detect_backend"]] = function(conn) {
    if (identical(conn@state[["server"]], "BigQuery"))
      dbreg:::backend_info_from_metadata(connection_class = "BigQueryConnection") else dbreg:::detect_backend(conn)
  }
  for (name in c("sql_model_matrix", "sql_model_matrix_internal", "get_column_info", "i_dictionaries")) {
    fn = get(name, envir = asNamespace("dbreg"))
    environment(fn) = scope
    scope[[name]] = fn
  }
  compile = scope[["sql_model_matrix"]]
  formulas = list(~i(text), ~i(key), ~i(adjacent), ~i(text, flag), ~i(day), ~0 + i(text))
  captures = list()
  for (server in c("Amazon Athena", "Presto", "Trino", "Unclassified DBMS")) {
    state[["server"]] = server
    for (index in seq_along(formulas)) {
      state[["queries"]] = character()
      formula = formulas[[index]]
      environment(formula) = baseenv()
      compiled = compile(formula, conn, "existing_backend_input",
        intercept = attr(stats::terms(formula), "intercept") == 1L)
      query = paste("SELECT", paste(paste(compiled[["select_exprs"]], "AS", DBI::dbQuoteIdentifier(conn, compiled[["aliases"]])), collapse = ","),
        "FROM existing_backend_input ORDER BY id")
      values = DBI::dbGetQuery(conn, query)
      captures[[paste(server, index, sep = "-")]] = list(backend = scope[["detect_backend"]](conn)[["name"]],
        compiled = compiled, queries = state[["queries"]], values = values)
    }
  }
  state[["server"]] = "BigQuery"
  state[["types"]] = c(text = "STRING", key = "INT64", adjacent = "FLOAT64", flag = "BOOL", day = "DATE")
  for (index in 1:4) {
    state[["queries"]] = character()
    formula = list(~i(text), ~0 + i(text), ~i(key), ~i(adjacent))[[index]]
    environment(formula) = baseenv()
    compiled = tryCatch(compile(formula, conn, "existing_backend_input",
      intercept = attr(stats::terms(formula), "intercept") == 1L), error = function(e) e)
    error = if (inherits(compiled, "error")) conditionMessage(compiled) else NULL
    values = if (is.null(error)) DBI::dbGetQuery(conn, paste("SELECT",
      paste(paste(compiled[["select_exprs"]], "AS", DBI::dbQuoteIdentifier(conn, compiled[["aliases"]])), collapse = ","),
      "FROM existing_backend_input ORDER BY id")) else NULL
    captures[[paste0("BigQuery-", index)]] = list(backend = scope[["detect_backend"]](conn)[["name"]],
      compiled = if (is.null(error)) compiled else NULL, queries = state[["queries"]], values = values, error = error)
  }
  removed = list()
  for (server in c("SQLite", "PostgreSQL")) {
    state[["server"]] = server
    state[["queries"]] = character()
    error = tryCatch({compile(~i(text), conn, "existing_backend_input"); NULL}, error = function(e) conditionMessage(e))
    removed[[server]] = list(error = error, queries = state[["queries"]])
  }
  list(data = data, captures = captures, removed = removed)
}
