source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_glm_wide.R", package = "dbreg"), local = TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbglm_test_wide(conn)
  expect_identical(dbglm_validation_tables(conn), character())
})
