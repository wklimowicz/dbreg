if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_sql_prediction_keys.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    dbreg_test_sql_prediction_keys(conn)
  })
}
