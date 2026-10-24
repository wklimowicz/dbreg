# Native engine validation is opt-in and uses the same cases and references.
if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_fe_edges.R", package = "dbreg"), local = TRUE)
  local({
    conn = DBI::dbConnect(odbc::databricks(),
      driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
      httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
    on.exit(DBI::dbDisconnect(conn))
    DBI::dbExecute(conn, "SET use_cached_result=false")
    dbreg_test_fe_edges(conn, "glm", verbose = TRUE)
    dbreg_test_fe_edges(conn, "ols", verbose = TRUE)
  })
}
