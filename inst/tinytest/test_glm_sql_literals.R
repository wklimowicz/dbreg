# Large FE lookup literals must survive dialect translation byte for byte.
# This size previously exhausted PCRE and silently truncated the covariance SQL.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  conn = dbglm_test_databricks_connection(native)
  payload = strrep("QUJDZGVmMDEyMzQ1Njc4OSsv", 300000L)
  literal = paste0("'", payload, "''FSUM isfinite EXCLUDE VARCHAR 東京'")
  quoted = paste0(literal, ",\"FSUM\"\"VARCHAR\",`EXCLUDE``FSUM`")
  sql = paste0("SELECT ", dbreg:::sql_sum(conn, "x"), ",CAST(x AS VARCHAR),", quoted, "\nFROM t")
  expected = paste0("SELECT SUM(x),CAST(x AS STRING),", quoted, "\nFROM t")
  expect_identical(dbreg:::backend_sql(conn, sql), expected)
  expect_identical(dbreg:::backend_sql(native, sql), sql)
  expect_identical(dbreg:::backend_sql(conn, paste0("SELECT ", literal)), paste0("SELECT ", literal))
  result = DBI::dbGetQuery(native, dbreg:::backend_sql(conn,
    paste0("SELECT LENGTH(", literal, ") AS n,", dbreg:::sql_sum(conn, "x"), " AS total FROM (VALUES (1),(2)) t(x)")))
  expect_equal(as.numeric(result[["n"]]), nchar(payload) + nchar("'FSUM isfinite EXCLUDE VARCHAR 東京"))
  expect_equal(as.numeric(result[["total"]]), 3)
  # Doubled quote/backtick escapes and empty literals remain distinct tokens.
  expect_identical(dbreg:::backend_sql(conn, "SELECT '', 'a''b', \"a\"\"b\", `a``b`, FSUM(x)"),
    "SELECT '', 'a''b', \"a\"\"b\", `a``b`, FSUM(x)")
})
