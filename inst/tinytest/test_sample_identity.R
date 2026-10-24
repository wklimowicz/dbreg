source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_sample_identity.R", package = "dbreg"), local = TRUE)
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2")), bigint = "integer64")
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_sample_identity(conn)
})
