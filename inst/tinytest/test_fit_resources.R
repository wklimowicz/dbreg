library(DBI)
local({
  # Delegate real SQL to DuckDB while injecting driver failures at DDL boundaries.
  methods::setClass("dbreg_resource_test", contains = "DBIConnection",
    slots = c(native = "ANY", state = "environment"))
  methods::setMethod("dbGetInfo", "dbreg_resource_test", function(dbObj, ...) list(dbms.name = "DuckDB"))
  methods::setMethod("dbIsValid", "dbreg_resource_test", function(dbObj, ...) DBI::dbIsValid(dbObj@native))
  methods::setMethod("dbQuoteIdentifier", c("dbreg_resource_test", "character"),
    function(conn, x, ...) DBI::dbQuoteIdentifier(conn@native, x, ...))
  methods::setMethod("dbListTables", "dbreg_resource_test", function(conn, ...) DBI::dbListTables(conn@native, ...))
  methods::setMethod("dbGetQuery", c("dbreg_resource_test", "character"),
    function(conn, statement, ...) DBI::dbGetQuery(conn@native, statement, ...))
  methods::setMethod("dbExecute", c("dbreg_resource_test", "character"), function(conn, statement, ...) {
    state = conn@state
    state[["calls"]] = c(state[["calls"]], statement)
    snapshot = grepl("^CREATE TEMP TABLE.*dbreg_sample_", statement)
    drop = grepl("^DROP TABLE.*dbreg_sample_", statement)
    if (drop && state[["fail_drop"]] > 0) {
      state[["fail_drop"]] = state[["fail_drop"]] - 1
      stop("injected drop failure")
    }
    if (grepl('^CREATE TEMP TABLE "export"', statement) && isTRUE(state[["race"]])) {
      DBI::dbExecute(conn@native, 'CREATE TEMP TABLE "export" AS SELECT 999 AS id')
    }
    value = DBI::dbExecute(conn@native, statement, ...)
    if (snapshot && isTRUE(state[["mutate"]])) DBI::dbExecute(conn@native, "DELETE FROM source")
    if (snapshot && isTRUE(state[["fail_create"]])) stop("injected post-create failure")
    value
  })
  methods::setMethod("dbDisconnect", "dbreg_resource_test", function(conn, ...) {
    conn@state[["disconnects"]] = conn@state[["disconnects"]] + 1L
    DBI::dbDisconnect(conn@native, ...)
  })
  native = DBI::dbConnect(duckdb::duckdb())
  on.exit({
    DBI::dbDisconnect(native, shutdown = TRUE)
    methods::removeClass("dbreg_resource_test")
  })
  state = new.env(parent = emptyenv())
  state[["calls"]] = character()
  state[["disconnects"]] = 0L
  state[["fail_drop"]] = 0
  state[["mutate"]] = state[["fail_create"]] = state[["race"]] = FALSE
  conn = methods::new("dbreg_resource_test", native = native, state = state)
  set.seed(11111)
  d = data.frame(id = 1:60, x = rnorm(60), w = runif(60, .2, 2))
  d[["y"]] = .6 * d[["x"]] + rnorm(60)
  reset_source = function() DBI::dbWriteTable(native, "source", d, overwrite = TRUE)
  fit = function(...) dbreg(y ~ x, conn = conn, table = "source", row_id = "id", sample_table = "export", ...)
  reset_source()
  # Mutating the caller source after snapshot creation cannot change the fit/export.
  state[["mutate"]] = TRUE
  model = fit()
  expect_equal(as.numeric(coef(model)), as.numeric(coef(lm(y ~ x, d))), tolerance = 1e-10)
  expect_equal(DBI::dbGetQuery(native, "SELECT COUNT(*) AS n FROM source")[["n"]], 0)
  expect_identical(DBI::dbGetQuery(native, "SELECT id FROM export ORDER BY id")[["id"]], d[["id"]])
  expect_equal(sort(DBI::dbListTables(native)), c("export", "source"))
  DBI::dbExecute(native, "DROP TABLE export")
  state[["mutate"]] = FALSE
  reset_source()
  # The driver may report an error after the private CREATE has actually succeeded.
  state[["fail_create"]] = TRUE
  expect_error(fit(), "post-create failure")
  expect_identical(DBI::dbListTables(native), "source")
  state[["fail_create"]] = FALSE
  DBI::dbExecute(native, "UPDATE source SET w=-1 WHERE id=1")
  expect_error(fit(weights = "w"), "non-negative")
  expect_identical(DBI::dbListTables(native), "source")
  reset_source()
  # A concurrent caller table wins: CREATE must fail without replacing or deleting it.
  state[["race"]] = TRUE
  expect_error(fit(), "already exists")
  expect_equal(DBI::dbGetQuery(native, "SELECT id FROM export")[["id"]], 999)
  expect_equal(sort(DBI::dbListTables(native)), c("export", "source"))
  DBI::dbExecute(native, "DROP TABLE export")
  state[["race"]] = FALSE
  # Failed cleanup remains tracked. Neither model removal nor GC calls the driver.
  state[["fail_drop"]] = 2
  model = fit()
  pending = DBI::dbListTables(native)
  expect_equal(sum(grepl("^dbreg_sample_", pending)), 1)
  calls = state[["calls"]]
  rm(model)
  invisible(gc())
  expect_identical(state[["calls"]], calls)
  expect_error(dbreg:::flush_retired_resources(conn), "retired fit resources")
  expect_identical(sort(DBI::dbListTables(native)), sort(pending))
  dbreg:::flush_retired_resources(conn)
  expect_equal(sort(DBI::dbListTables(native)), c("export", "source"))
  expect_identical(DBI::dbGetQuery(native, "SELECT id FROM export ORDER BY id")[["id"]], d[["id"]])
  DBI::dbExecute(native, "DROP TABLE export")
  state[["fail_drop"]] = 1
  model = fit()
  DBI::dbExecute(native, "DROP TABLE export")
  # The next fit retries synchronously, before allocating another snapshot.
  model = fit()
  expect_equal(sort(DBI::dbListTables(native)), c("export", "source"))
  expect_identical(state[["disconnects"]], 0L)
  expect_true(DBI::dbIsValid(native))
  expect_equal(length(ls(dbreg:::.dbreg_retired_resources, all.names = TRUE)), 0)
  # Owned connections close once; a second release is harmless.
  owned = DBI::dbConnect(duckdb::duckdb())
  resources = dbreg:::fit_resources()
  resources[["conn"]] = owned
  resources[["own_conn"]] = TRUE
  dbreg:::release_fit_resources(resources, quiet = FALSE)
  expect_false(DBI::dbIsValid(owned))
  dbreg:::release_fit_resources(resources, quiet = FALSE)
  expect_false(dbreg:::resource_pending(resources))
})
