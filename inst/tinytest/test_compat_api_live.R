if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")),"true")) {
  source(system.file("tinytest/helper_glm_reference.R",package="dbreg"),local=TRUE)
  source(system.file("tinytest/helper_glm_validation.R",package="dbreg"),local=TRUE)
  source(system.file("tinytest/helper_compat_api.R",package="dbreg"),local=TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    dbreg_test_compat_api(conn)
  })
}
