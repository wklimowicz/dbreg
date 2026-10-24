source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_sql_prediction_keys.R", package = "dbreg"), local = TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_sql_prediction_keys(conn)
})
