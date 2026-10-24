if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_linear_design.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(live = TRUE)
    on.exit(DBI::dbDisconnect(conn))
    dbreg_test_linear_backend(conn)
  })
}
