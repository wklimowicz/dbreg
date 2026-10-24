source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_linear_design.R", package = "dbreg"), local = TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_linear_backend(conn)
})
