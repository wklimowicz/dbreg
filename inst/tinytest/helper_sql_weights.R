dbreg_test_sql_weights = function(conn) {
  before = dbglm_validation_tables(conn)
  tables = character()
  on.exit(for (table in rev(tables)) dbglm_validation_drop(conn, table))
  for (value in c("Infinity", "-Infinity")) {
    table = basename(tempfile("glm_validation_weights_"))
    tables = c(tables, table)
    DBI::dbExecute(conn, paste0("CREATE TEMP TABLE ", table,
      " AS SELECT 1.0 AS y,2.0 AS x,3.0 AS endog,4.0 AS instrument,CAST('", value, "' AS DOUBLE) AS w"))
    expect_error(dbreg(y ~ x, conn = conn, table = table, weights = "w"), "Weights must be (finite|non-negative)")
    expect_error(dbivreg(y ~ x | endog ~ instrument, conn = conn, table = table, weights = "w"), "Weights must be (finite|non-negative)")
    expect_true(DBI::dbIsValid(conn))
    expect_identical(dbglm_validation_tables(conn), sort(c(before, tables)))
  }
  for (table in tables) dbglm_validation_drop(conn, table)
  expect_identical(dbglm_validation_tables(conn), before)
}
