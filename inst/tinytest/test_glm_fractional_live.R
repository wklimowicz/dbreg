if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_fractional.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    dbglm_test_fractional_weights(conn, verbose = TRUE)
    dbglm_test_fractional_scale(conn, verbose = TRUE)
    dbglm_test_fractional_mutation(conn, verbose = TRUE)
  })
}
