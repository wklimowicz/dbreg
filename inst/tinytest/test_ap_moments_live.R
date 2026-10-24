if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_ap_moments.R", package = "dbreg"), local = TRUE)
  local({
    conn = DBI::dbConnect(odbc::databricks(),
      driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
      httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
    on.exit(DBI::dbDisconnect(conn))
    DBI::dbExecute(conn, "SET use_cached_result=false")
    dbreg_test_ap_moments(conn, verbose = TRUE)
  })
}
