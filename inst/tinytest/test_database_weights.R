library(DBI)

# Database sources must reject infinite weights before estimation, just as
# in-memory inputs do. Missing/zero weights retain their existing semantics.
local({
  set.seed(721)
  n = 80L
  d = data.frame(x = rnorm(n), z = rnorm(n), w = rep(1, n))
  d[["endog"]] = .6*d[["z"]] + .3*d[["x"]] + rnorm(n)
  d[["y"]] = 1 + .5*d[["x"]] + .4*d[["endog"]] + rnorm(n)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "weight_source", d)
  inventory = sort(DBI::dbListTables(conn))
  fits = list(
    function() dbreg(y ~ x, conn = conn, table = "weight_source", weights = "w", strategy = "moments"),
    function() dbivreg(y ~ x | endog ~ z, conn = conn, table = "weight_source", weights = "w", strategy = "moments")
  )
  for (value in c(Inf, -Inf)) {
    bad = d
    bad[["w"]][1L] = value
    DBI::dbWriteTable(conn, "weight_source", bad, overwrite = TRUE)
    expect_error(dbreg(y ~ x, data = bad, weights = "w"), "Weights must be finite")
    for (fit in fits) {
      expect_error(fit(), "Weights must be (finite|non-negative)")
      expect_true(DBI::dbIsValid(conn))
      expect_identical(sort(DBI::dbListTables(conn)), inventory)
    }
  }
  d[["w"]][1:2] = c(0, NA_real_)
  DBI::dbWriteTable(conn, "weight_source", d, overwrite = TRUE)
  expected = stats::lm(y ~ x, data = d, subset = w > 0)
  actual = fits[[1L]]()
  expect_equal(unname(coef(actual)), unname(coef(expected)), tolerance = 1e-10)
  expect_equal(as.numeric(vcov(actual)), as.numeric(vcov(expected)), tolerance = 1e-10)
  expect_identical(sort(DBI::dbListTables(conn)), inventory)
  d[["w"]][1L] = -.1
  DBI::dbWriteTable(conn, "weight_source", d, overwrite = TRUE)
  for (fit in fits) expect_error(fit(), "Weights must be non-negative")
})

# Exercise Spark's generated validation predicate through a local DBI fixture.
# This checks SQL semantics, not Spark execution/ODBC; a live probe is separate.
local({
  if (!methods::isClass("dbreg_weight_sql_fixture")) {
    methods::setClass("dbreg_weight_sql_fixture", contains = "DBIConnection", slots = c(native = "ANY"))
    methods::setMethod("dbGetInfo", "dbreg_weight_sql_fixture", function(dbObj, ...) list(dbms.name = "spark"))
    methods::setMethod("dbGetQuery", c("dbreg_weight_sql_fixture", "character"),
      function(conn, statement, ...) DBI::dbGetQuery(methods::slot(conn, "native"), statement, ...))
  }
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  spark = methods::new("dbreg_weight_sql_fixture", native = native)
  for (value in c("Infinity", "-Infinity")) {
    from = paste0("FROM (SELECT CAST('", value, "' AS DOUBLE) AS w) source")
    expect_error(dbreg:::validate_weights("w", NULL, spark, from), "Weights must be (finite|non-negative)")
  }
  # SQL NaN remains governed by existing missing-value behavior; it is not
  # reclassified as an infinite weight by this targeted validation change.
  for (value in c("1.0", "0.0", "NULL", "CAST('NaN' AS DOUBLE)")) {
    expect_silent(dbreg:::validate_weights("w", NULL, spark, paste0("FROM (SELECT ", value, " AS w) source")))
  }
  expect_error(dbreg:::validate_weights("w", NULL, spark, "FROM (SELECT -1 AS w) source"), "Weights must be non-negative")
})
