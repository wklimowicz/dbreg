# Opt-in: credentials/driver must already be configured. The same numerical
# cases run offline too; this gate verifies the actual Spark/ODBC execution.
if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  if (!requireNamespace("odbc", quietly = TRUE)) stop("Live Databricks tests require odbc.")
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_databricks_cases.R", package = "dbreg"), local = TRUE)
  local({
    conn = DBI::dbConnect(odbc::databricks(),
      driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
      httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
    on.exit(DBI::dbDisconnect(conn))
    DBI::dbExecute(conn, "SET use_cached_result=false")
    dbglm_test_databricks_cases(conn, verbose = TRUE)
  })
}
