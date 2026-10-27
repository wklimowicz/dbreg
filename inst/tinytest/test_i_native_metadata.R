library(DBI)
# The fixture supplies native result metadata and records the original dialect
# SQL. DuckDB evaluates the resulting indicator values; this is not a live
# SQL Server/MySQL driver certification.
local({
  old = options(dbreg.i_max_columns = 1000L)
  on.exit(options(old))
  connection_class = "dbreg_i_metadata_connection"
  result_class = "dbreg_i_metadata_result"
  if (!methods::isClass(connection_class)) methods::setClass(connection_class,
    contains = "DBIConnection", slots = c(native = "ANY", state = "environment"))
  if (!methods::isClass(result_class)) methods::setClass(result_class,
    contains = "DBIResult", slots = c(native = "ANY", state = "environment"))
  translate = function(statement) {
    limit = regmatches(statement, regexec("^SELECT TOP ([0-9]+) ", statement))[[1L]]
    if (length(limit)) statement = paste(sub("^SELECT TOP [0-9]+ ", "SELECT ", statement), "LIMIT", limit[2L])
    statement = gsub("VARCHAR(MAX)", "VARCHAR", statement, fixed = TRUE)
    statement = gsub(" AS SIGNED)", " AS BIGINT)", statement, fixed = TRUE)
    statement = gsub("`([^`]+)`", '"\\1"', statement)
    gsub("\\[([^]]+)\\]", '"\\1"', statement)
  }
  record = function(conn, statement) {
    if (isTRUE(conn@state[["active"]])) stop("A schema result is still open.")
    conn@state[["sql"]] = c(conn@state[["sql"]], statement)
    if (grepl("typeof(", statement, fixed = TRUE)) stop("Native driver fixture rejects typeof().")
    if (conn@state[["server"]] == "Microsoft SQL Server" && grepl(" LIMIT ", statement, fixed = TRUE)) {
      stop("Native SQL Server fixture rejects LIMIT.")
    }
  }
  methods::setMethod("dbGetInfo", connection_class, function(dbObj, ...) {
    list(dbms.name = dbObj@state[["server"]])
  })
  methods::setMethod("dbQuoteIdentifier", c(connection_class, "character"), function(conn, x, ...) {
    if (conn@state[["server"]] == "Microsoft SQL Server") DBI::SQL(paste0("[", gsub("]", "]]", x, fixed = TRUE), "]")) else
      DBI::SQL(paste0("`", gsub("`", "``", x, fixed = TRUE), "`"))
  })
  methods::setMethod("dbQuoteLiteral", c(connection_class, "character"), function(conn, x, ...) {
    DBI::dbQuoteLiteral(conn@native, x, ...)
  })
  methods::setMethod("dbSendQuery", c(connection_class, "character"), function(conn, statement, ...) {
    record(conn, statement)
    native = DBI::dbSendQuery(conn@native, translate(statement), ...)
    conn@state[["active"]] = TRUE
    methods::new(result_class, native = native, state = conn@state)
  })
  methods::setMethod("dbColumnInfo", result_class, function(res, ...) {
    if (isTRUE(res@state[["fail_info"]])) stop("Column metadata failure.")
    res@state[["columns"]]
  })
  methods::setMethod("dbFetch", result_class, function(res, n = -1, ...) {
    if (isTRUE(res@state[["fail_fetch"]])) stop("Schema fetch failure.")
    DBI::dbFetch(res@native, n = n, ...)
  })
  methods::setMethod("dbClearResult", result_class, function(res, ...) {
    res@state[["active"]] = FALSE
    res@state[["cleared"]] = res@state[["cleared"]] + 1L
    DBI::dbClearResult(res@native, ...)
  })
  methods::setMethod("dbGetQuery", c(connection_class, "character"), function(conn, statement, ...) {
    record(conn, statement)
    if (grepl("SQL_VARIANT_PROPERTY", statement, fixed = TRUE)) return(conn@state[["probe"]])
    if (grepl("AS native_decimal_text", statement, fixed = TRUE)) return(conn@state[["probe"]])
    if (grepl(" GROUP BY ", statement, fixed = TRUE) && !is.null(conn@state[["keys"]])) {
      return(data.frame(value = conn@state[["keys"]]))
    }
    DBI::dbGetQuery(conn@native, translate(statement), ...)
  })
  methods::setMethod("dbExecute", c(connection_class, "character"), function(conn, statement, ...) {
    stop("Metadata/design compilation must not execute physical DDL.")
  })
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  state = new.env(parent = emptyenv())
  state[["server"]] = "Microsoft SQL Server"
  state[["columns"]] = data.frame(name = c("key", "flag", "x"), type = c("-5", "-7", "8"))
  state[["sql"]] = character()
  state[["cleared"]] = 0L
  state[["active"]] = FALSE
  conn = methods::new(connection_class, native = native, state = state)
  data = data.frame(key = rep(c("9007199254740992", "9007199254740993", "9007199254740994"), 4),
    flag = rep(c(FALSE, TRUE), 6), x = seq_len(12) / 7)
  DBI::dbWriteTable(native, "i_data_text", data)
  DBI::dbExecute(native, "CREATE TABLE i_data AS SELECT CAST(key AS BIGINT) AS key, flag, x FROM i_data_text")

  # Result cleanup is exercised through real DBI result methods, including both
  # failure boundaries. No schema result survives into a dictionary SELECT.
  query = "SELECT TOP 0 key, flag, x FROM i_data"
  schema = dbreg:::read_i_schema(conn, query)
  expect_identical(state[["cleared"]], 1L)
  expect_false(state[["active"]])
  for (failure in c("fail_info", "fail_fetch")) {
    state[[failure]] = TRUE
    expect_error(dbreg:::read_i_schema(conn, query), "failure")
    expect_false(state[["active"]])
    state[[failure]] = FALSE
  }
  expect_identical(state[["cleared"]], 3L)

  # Only the driver's class identification is simulated. Metadata contents,
  # scoped schema fetch, dictionary SQL and value evaluation cross actual DBI.
  scope = new.env(parent = asNamespace("dbreg"))
  reader = function(conn, query) {
    result = dbreg:::read_i_schema(conn, query)
    result[["driver"]] = "odbc"
    result
  }
  scope[["read_i_schema"]] = reader
  get_info = dbreg:::get_column_info
  environment(get_info) = scope
  scope[["get_column_info"]] = get_info
  # The public compiler delegates through this owner on the certified prefix.
  internal = dbreg:::sql_model_matrix_internal
  environment(internal) = scope
  scope[["sql_model_matrix_internal"]] = internal
  compile = dbreg::sql_model_matrix
  environment(compile) = scope
  declaration = dbreg:::parse_sql_i_term(quote(i(key, flag)), baseenv())
  info = get_info(conn, "i_data", c("key", "flag"), level_vars = "key", i_declarations = list(declaration))
  expect_identical(info[["classes"]][["flag"]], "logical")
  expect_identical(info[["logical_sql"]][["flag"]], "[flag] = 1")
  expect_identical(info[["i_specs"]][[1L]][["values"]], unique(data[["key"]]))
  expect_identical(dbreg:::sql_value_literals(conn, c(TRUE, FALSE, NA)), c("1", "0", "NULL"))
  compiled = compile(~i(key, flag), conn, "i_data")
  expect_true(all(grepl("AND ([flag] = 1)", compiled[["select_exprs"]], fixed = TRUE)))
  actual = DBI::dbGetQuery(conn, paste("SELECT", paste(compiled[["select_exprs"]], collapse = ","), "FROM i_data"))
  expected = outer(data[["key"]], unique(data[["key"]]), `==`) * data[["flag"]]
  expect_equal(unname(as.matrix(actual)), unname(expected))
  expect_equal(unname(dbreg:::evaluate_i_spec(info[["i_specs"]][[1L]], data[["key"]], data[["flag"]])), unname(expected))
  expect_true(any(grepl("SELECT TOP 1002 CAST([key] AS VARCHAR(MAX))", state[["sql"]], fixed = TRUE)))
  expect_false(any(grepl("LIMIT", state[["sql"]], fixed = TRUE)))

  # SQL prediction recompilation must retain BIT predicates and frozen levels.
  object = list(model_spec = list(rhs = ~i(key, flag), source_vars = c("key", "flag"),
    factor_levels = info[["levels"]], factor_classes = info[["factor_classes"]],
    i_specs = info[["i_specs"]], source_classes = info[["classes"]],
    columns = data.frame(name = compiled[["coef_names"]])), fe = character(), has_intercept = TRUE)
  class(object) = "dbglm"
  design = dbreg:::sql_fitted_prediction_design(object, conn)
  expect_true(all(grepl("AND ([flag] = 1)", design[["expressions"]][-1L], fixed = TRUE)))

  # Decimal/FLOAT probes use native evidence only, over a quoted filtered source.
  schema[["driver"]] = "odbc"
  schema[["columns"]] = data.frame(name = "key", type = "3")
  state[["probe"]] = data.frame(base_type = "decimal", precision = 38L, scale = 0L)
  metadata = dbreg:::i_native_metadata(conn, "(SELECT key FROM i_data WHERE x > 0) src", "key", schema)
  expect_identical(metadata, list(types = "integer", sql_types = "DECIMAL(38,0)"))
  expect_true(grepl("SELECT TOP 1 CAST(SQL_VARIANT_PROPERTY", tail(state[["sql"]], 1), fixed = TRUE))
  expect_true(grepl("WHERE [key] IS NOT NULL", tail(state[["sql"]], 1), fixed = TRUE))
  state[["probe"]][["scale"]] = 2L
  expect_error(dbreg:::i_native_metadata(conn, "#source", "key", schema), "fractional DECIMAL")
  schema[["columns"]][["type"]] = "6"
  for (precision in c(24L, 53L)) {
    state[["probe"]] = data.frame(base_type = "float", precision = precision, scale = 0L)
    expect_identical(dbreg:::i_native_metadata(conn, "#source", "key", schema)[["types"]],
      if (precision == 24L) "float" else "numeric")
  }
  state[["probe"]][["precision"]] = NA_integer_
  expect_error(dbreg:::i_native_metadata(conn, "#source", "key", schema), "invalid native numeric")
  for (code in c("12", "-9")) {
    schema[["columns"]][["type"]] = code
    before = length(state[["sql"]])
    expect_identical(dbreg:::i_native_metadata(conn, "#source", "key", schema)[["types"]], "character")
    expect_identical(length(state[["sql"]]), before)
  }
  schema[["columns"]][["type"]] = "-150"
  expect_error(dbreg:::i_native_metadata(conn, "#source", "key", schema), "unsupported native metadata")

  state[["server"]] = "MySQL"
  schema[["columns"]][["type"]] = "3"
  for (text in c("0", "0000000000001", "-9223372036854775808", "9007199254740993")) {
    state[["probe"]] = data.frame(native_decimal_text = text)
    expect_identical(dbreg:::i_native_metadata(conn, "i_data", "key", schema),
      list(types = "integer", sql_types = "DECIMAL"))
  }
  for (text in c("0.00", "1.00", "-1.20")) {
    state[["probe"]] = data.frame(native_decimal_text = text)
    expect_error(dbreg:::i_native_metadata(conn, "i_data", "key", schema), "fractional DECIMAL")
  }
  for (text in list("1e3", "1.", " 1", NA_character_, character(), 1)) {
    state[["probe"]] = data.frame(native_decimal_text = text)
    expect_error(dbreg:::i_native_metadata(conn, "i_data", "key", schema), "native DECIMAL text scale")
  }
  expect_true(grepl("LIMIT 1$", tail(state[["sql"]], 1)))
  expect_false(any(grepl("typeof", state[["sql"]], fixed = TRUE)))
  literal = dbreg:::sql_factor_literals(conn, list(key = "9007199254740993"), c(key = "integer64"))
  expect_identical(literal[["key"]], "CAST('9007199254740993' AS SIGNED)")
  state[["probe"]] = data.frame(native_decimal_text = "0")
  state[["columns"]] = data.frame(name = c("key", "flag"), type = c("3", "-7"))
  decimal = get_info(conn, "i_data", "key", level_vars = "key", i_declarations = list(dbreg:::parse_sql_i_term(quote(i(key)), baseenv())))
  expect_identical(decimal[["i_specs"]][[1L]][["values"]], unique(data[["key"]]))
  expect_equal(unname(dbreg:::evaluate_i_spec(decimal[["i_specs"]][[1L]], data[["key"]])),
    unname(outer(data[["key"]], unique(data[["key"]])[-1L], `==`) * 1))

  for (keys in list(c("000000000000000000001", "0002"),
                    c("-9223372036854775808", "9223372036854775807"))) {
    state[["keys"]] = keys
    dictionary = get_info(conn, "i_data", "key", level_vars = "key",
      i_declarations = list(dbreg:::parse_sql_i_term(quote(i(key)), baseenv())))
    expected_keys = sub("^[0-9]+:", "", dbreg:::prediction_key_part(keys, "integer"))
    expect_identical(dictionary[["i_specs"]][[1L]][["values"]], expected_keys)
  }
  for (keys in list("9223372036854775808", "-9223372036854775809", 9007199254740992)) {
    state[["keys"]] = keys
    expect_error(get_info(conn, "i_data", "key", level_vars = "key",
      i_declarations = list(dbreg:::parse_sql_i_term(quote(i(key)), baseenv()))),
      if (is.character(keys)) "signed BIGINT" else "native character transport")
  }
  state[["keys"]] = NULL
  # Ordinary numeric OLS/IV schema discovery still uses the original DBI route.
  before = state[["cleared"]]
  plain = get_info(conn, "i_data", "x")
  expect_identical(plain[["classes"]], c(x = "numeric"))
  expect_identical(state[["cleared"]], before)

  # Generic character keys retain upstream server ordering, even when an ENUM
  # marker is unavailable. Neither default reference nor keep re-sorts levels.
  DBI::dbExecute(native, "CREATE TYPE native_enum AS ENUM ('z', 'a', 'm')")
  DBI::dbExecute(native, "CREATE TABLE enum_data AS SELECT CAST(v AS native_enum) AS key FROM (VALUES ('m'),('z'),('a')) s(v)")
  state[["columns"]] = data.frame(name = "key", type = "12")
  for (call in list(quote(i(key)), quote(i(key, keep = c("m", "a"))))) {
    decl = dbreg:::parse_sql_i_term(call, baseenv())
    enum = get_info(conn, "enum_data", "key", level_vars = "key", i_declarations = list(decl))
    spec = enum[["i_specs"]][[1L]]
    expect_identical(as.character(spec[["values"]]), c("z", "a", "m"))
    expect_identical(spec[["names"]], c("key::a", "key::m"))
  }
  schema[["driver"]] = "rmysql"
  for (label in c("BIGINT", "DOUBLE", "STRING")) {
    schema[["columns"]][["type"]] = label
    expect_identical(dbreg:::i_native_metadata(conn, "i_data", "key", schema)[["types"]],
      switch(label, BIGINT = "integer", DOUBLE = "numeric", STRING = "character"))
  }
  schema[["columns"]][["type"]] = "<unknown>"
  expect_error(dbreg:::i_native_metadata(conn, "i_data", "key", schema), "unsupported native metadata")
  schema[["driver"]] = "odbc"
  schema[["columns"]][["type"]] = "3"
  state[["server"]] = "MariaDB"
  expect_error(dbreg:::i_native_metadata(conn, "i_data", "key", schema), "unverified")
  expect_false(state[["active"]])
})
