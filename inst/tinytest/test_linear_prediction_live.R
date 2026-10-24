if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_linear_prediction.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    dbreg_test_linear_prediction(conn)
    dbreg_test_linear_prediction_keys(conn)
    dbreg_test_linear_prediction_estimability(conn)
    dbreg_test_cluster_key_precision(conn)
    dbreg_test_linear_prediction_snapshot(conn)
  })
}
