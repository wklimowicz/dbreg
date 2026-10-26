# Construction owns finite predicates; remaining dialect adapter does not.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", temp_directory = "")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  remote = dbglm_test_databricks_connection(native)
  expect_identical(dbreg:::sql_is_finite(native, c("x", "(x+y)")), c("isfinite(x)", "isfinite((x+y))"))
  expect_identical(dbreg:::sql_is_finite(native, character()), character())
  # Golden outputs captured from the former Spark finite-token expansion.
  expected = c(
    "(NOT isnan(CAST((x) AS DOUBLE)) AND ABS(CAST((x) AS DOUBLE)) <= 1.7976931348623157e308)",
    "(NOT isnan(CAST(((x+y)) AS DOUBLE)) AND ABS(CAST(((x+y)) AS DOUBLE)) <= 1.7976931348623157e308)"
  )
  expect_identical(dbreg:::sql_is_finite(remote, c("x", "(x+y)")), expected)
  # The stabilized sum's finite guard changes ownership, not SUM/DECIMAL policy.
  legacy_guard = paste0(
    "(NOT isnan(CAST(((x)) AS DOUBLE)) AND ABS(CAST(((x)) AS DOUBLE)) <= 1.7976931348623157e308)",
    " AND ABS((x)) < 1e10"
  )
  legacy_sum = paste0(
    "(CASE WHEN MAX(ABS(x)) < 1e-9 THEN SUM(x) WHEN SUM(CASE WHEN ", legacy_guard,
    " THEN 0 ELSE 1 END) = 0 THEN CAST(SUM(CAST(CASE WHEN ", legacy_guard,
    " THEN (x) END AS DECIMAL(38,18))) AS DOUBLE) ELSE SUM(x) END)"
  )
  expect_identical(dbreg:::sql_stable_sum(remote, "x"), legacy_sum)
  expect_identical(dbreg:::backend_sql(remote, dbreg:::sql_stable_sum(remote, "x")), legacy_sum)
  expect_identical(dbreg:::backend_sql(remote, "SELECT isfinite(x),'isfinite(x)' AS label"),
    "SELECT isfinite(x),'isfinite(x)' AS label")
  expect_error(dbreg:::sql_is_finite(NULL, "x", list(name = "mysql")), "DuckDB and Spark only")
  # Independent SQL truth table: NULL remains unknown; NaN and infinities false.
  fixture = paste0("(SELECT id,v,0.0 AS o FROM (VALUES ",
    "(1,CAST(NULL AS DOUBLE)),(2,CAST('NaN' AS DOUBLE)),",
    "(3,CAST('Infinity' AS DOUBLE)),(4,CAST('-Infinity' AS DOUBLE)),",
    "(5,0.0),(6,-0.0),(7,1e300),(8,1e-300)) t(id,v))")
  expected_finite = c(NA, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  for (conn in list(native, remote)) {
    actual = dbreg:::backend_get_query(conn, paste("SELECT", dbreg:::sql_is_finite(conn, "v"),
      "AS finite FROM", fixture, "ORDER BY id"))
    expect_identical(actual[["finite"]], expected_finite)
    # The real offset builder must preserve the same ordered truth values.
    offset = dbreg:::offset_source_frame(conn, paste("FROM", fixture), "v", "o")
    actual = dbreg:::backend_get_query(conn, paste("SELECT", offset[["finite"]],
      "AS finite FROM", fixture, "ORDER BY id"))
    expect_identical(actual[["finite"]], expected_finite)
    # Validate original weights before any outcome-missing filter; NaN itself
    # remains missing, while an infinite weight on a missing-outcome row errors.
    expect_error(dbreg:::validate_weights("w", NULL, conn,
      "FROM (SELECT CAST('Infinity' AS DOUBLE) AS w,CAST(NULL AS DOUBLE) AS y) t"), "Weights must be finite")
    expect_error(dbreg:::validate_weights("w", NULL, conn,
      "FROM (SELECT CAST('-Infinity' AS DOUBLE) AS w) t"), "Weights must be non-negative")
    expect_silent(dbreg:::validate_weights("w", NULL, conn,
      "FROM (SELECT CAST('NaN' AS DOUBLE) AS w UNION ALL SELECT CAST(NULL AS DOUBLE)) t"))
    # Real absorbed-GLM diagnostics count only the three invalid finite eta
    # rows. NULL eta is already excluded upstream and retains SQL unknown here.
    working = paste0("(SELECT *,1.0 AS w,1.0 AS x,0.0 AS z,0.0 AS score,0.0 AS dev,",
      "0.5 AS mu,v AS eta FROM ", fixture, ")")
    diagnostic = dbreg:::dbglm_fe_evaluate(conn, working, "x", character())
    expect_equal(as.numeric(diagnostic[["bad"]]), 3)
    # Existing SQL-row projection certificate must reject invalid rows rather
    # than accepting a zero/NaN mean as a successful projection.
    invalid = paste0("(SELECT 1 AS g,1.0 AS w,v FROM ", fixture, ")")
    expect_identical(dbreg:::dbreg_projection_certificate(conn, invalid, "g", "v", 8), Inf)
    valid = "(SELECT 1 AS g,1.0 AS w,v FROM (VALUES (-1.0),(1.0)) t(v))"
    expect_equal(dbreg:::dbreg_projection_certificate(conn, valid, "g", "v", 2), 0)
  }
})
