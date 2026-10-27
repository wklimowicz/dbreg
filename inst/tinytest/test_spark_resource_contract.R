library(DBI)
# Capability checks dominate physical resource acquisition and AP cleanup.
local({
  state = new.env(parent = emptyenv())
  state[["events"]] = character()
  record = function(name) function(...) {
    state[["events"]] = c(state[["events"]], name)
    stop("Resource dispatch must not occur.", call. = FALSE)
  }
  project = dbreg:::dbreg_alternating_projections
  environment(project) = new.env(parent = environment(project))
  environment(project)[["detect_backend"]] = function(conn) list(name = "spark",
    is_databricks = TRUE, supports_physical_temp_tables = TRUE, native_query = FALSE)
  for (name in c("create_temp_table_as", "drop_table_if_exists",
                 "dbreg_moment_projections", "dbreg_varying_projections")) {
    environment(project)[[name]] = record(name)
  }
  old = options(dbreg.ap_solver = "iterative", dbreg.ap_method = "join")
  on.exit(options(old))
  inputs = list2env(list(conn = NULL, from_statement = "FROM source", fe = c("g", "h"),
    yvar = "y", weights = NULL, cluster_var = NULL, verbose = FALSE))
  for (method in c("join", "window")) {
    options(dbreg.ap_method = method)
    expect_error(project(inputs, "x", "x", 100L, 1e-10), "native query mode")
    expect_identical(state[["events"]], character())
  }
  inputs[["fe_slopes"]] = list(g = "x")
  expect_error(project(inputs, "x", "x", 100L, 1e-10), "native query mode")
  expect_identical(state[["events"]], character())
  expect_error(project(inputs, "x", "x", 0, 1e-10), "ap_max_iter")
  expect_error(project(inputs, "x", "x", 100L, 0), "ap_tol")
  expect_identical(state[["events"]], character())
})

# A DBI boundary fixture delegates accepted DDL to DuckDB. Only the fixture's
# DROP TEMP/managed-Parquet syntax is adapted; original statements are recorded.
local({
  class_name = "dbreg_spark_resource_fixture"
  if (!methods::isClass(class_name)) methods::setClass(class_name,
    contains = "DBIConnection", slots = c(native = "ANY", state = "environment", metadata = "list"))
  methods::setMethod("dbGetInfo", class_name, function(dbObj, ...) dbObj@metadata)
  methods::setMethod("dbIsValid", class_name, function(dbObj, ...) DBI::dbIsValid(dbObj@native))
  methods::setMethod("dbQuoteIdentifier", c(class_name, "character"), function(conn, x, ...) {
    DBI::dbQuoteIdentifier(conn@native, x, ...)
  })
  methods::setMethod("dbGetQuery", c(class_name, "character"), function(conn, statement, ...) {
    DBI::dbGetQuery(conn@native, statement, ...)
  })
  methods::setMethod("dbExecute", c(class_name, "character"), function(conn, statement, ...) {
    conn@state[["statements"]] = c(conn@state[["statements"]], statement)
    if (isTRUE(conn@state[["forbid_execute"]])) stop("Physical DBI execution was forbidden.", call. = FALSE)
    statement = sub("^DROP TEMP TABLE", "DROP TABLE", statement)
    statement = gsub(" USING PARQUET", "", statement, fixed = TRUE)
    DBI::dbExecute(conn@native, statement, ...)
  })
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  state = new.env(parent = emptyenv())
  state[["statements"]] = character()
  state[["forbid_execute"]] = TRUE
  conn = methods::new(class_name, native = native, state = state,
    metadata = list(dbms.name = "Databricks SQL"))
  attr(conn, "connectionString") = list(useNativeQuery = FALSE)
  expect_error(dbreg:::create_temp_table_as(conn, "guard_invalid", "SELECT 1 AS x", "spark"), "native query mode")
  expect_identical(state[["statements"]], character())
  expect_identical(DBI::dbListTables(native), character())

  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  manager = dbreg:::fit_table_manager(conn, resources)
  old = options(dbreg.glm_databricks_materialize_intermediates = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(manager[["create"]]("SELECT 1 AS x", materialize = FALSE), "(SELECT 1 AS x)")
  expect_error(manager[["create"]]("SELECT 1 AS x"), "native query mode")
  expect_error(manager[["create"]]("SELECT 1 AS x", target = "existing_target"), "native query mode")
  expect_identical(resources[["tables"]], character())
  expect_identical(resources[["views"]], character())
  dbreg:::release_fit_resources(resources)
  expect_identical(state[["statements"]], character())
  expect_true(DBI::dbIsValid(native))

  # TRUE and unknown native metadata preserve the accepted session-table path.
  state[["forbid_execute"]] = FALSE
  for (native_mode in list(TRUE, NULL)) {
    attr(conn, "connectionString") = if (is.null(native_mode)) NULL else list(useNativeQuery = native_mode)
    dbreg:::create_temp_table_as(conn, "guard_valid", "SELECT 7 AS x", "spark")
    expect_equal(DBI::dbGetQuery(native, "SELECT x FROM guard_valid")[["x"]], 7)
    expect_true(startsWith(utils::tail(state[["statements"]], 1L), "CREATE OR REPLACE TEMPORARY TABLE"))
    dbreg:::drop_table_if_exists(conn, "guard_valid", "spark")
    expect_identical(DBI::dbListTables(native), character())
  }
  oss = methods::new(class_name, native = native, state = state,
    metadata = list(dbms.name = "Apache Spark SQL"))
  dbreg:::create_temp_table_as(oss, "guard_managed", "SELECT 9 AS x", "spark")
  expect_equal(DBI::dbGetQuery(native, "SELECT x FROM guard_managed")[["x"]], 9)
  expect_true(grepl("CREATE TABLE guard_managed USING PARQUET", utils::tail(state[["statements"]], 1L), fixed = TRUE))
  dbreg:::drop_table_if_exists(oss, "guard_managed", "spark")
  expect_identical(DBI::dbListTables(native), character())

  # Refusing physical AP resources does not refuse read-only ordinary fitting.
  data = data.frame(x = seq_len(30) / 30)
  data[["y"]] = 1 + 2 * data[["x"]] + sin(seq_len(30)) / 10
  DBI::dbWriteTable(native, "guard_source", data)
  attr(conn, "connectionString") = list(useNativeQuery = FALSE)
  state[["statements"]] = character()
  state[["forbid_execute"]] = TRUE
  fit = dbreg::dbreg(y ~ x, conn = conn, table = "guard_source", strategy = "moments")
  reference = stats::lm(y ~ x, data = data)
  expect_equal(as.numeric(stats::coef(fit)), as.numeric(stats::coef(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-8)
  close(fit)
  expect_identical(state[["statements"]], character())
  expect_identical(DBI::dbListTables(native), "guard_source")
  expect_true(DBI::dbIsValid(native))
  DBI::dbRemoveTable(native, "guard_source")
})
