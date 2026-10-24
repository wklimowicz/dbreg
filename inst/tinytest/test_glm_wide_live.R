if (identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")) {
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
  source(system.file("tinytest/helper_glm_wide.R", package = "dbreg"), local = TRUE)
  local({
    conn = dbglm_validation_connection(TRUE)
    on.exit(DBI::dbDisconnect(conn))
    before = dbglm_validation_tables(conn)
    dbglm_test_wide(conn, verbose = TRUE)
    expect_identical(dbglm_validation_tables(conn), before)
  })
}
