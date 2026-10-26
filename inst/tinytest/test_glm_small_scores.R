# Exercise the production Spark score expression through the local SQL fixture.
# Native learned-effect tests also check its effect on actual fitted predictors.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  native = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  remote = dbglm_test_databricks_connection(native)
  for (values in c("(1e-20),(2e-20)", "(1e-100),(2e-100)", "(1e-20),(-1e-20),(1e-35)", "(1.0),(2.0)")) {
    row = dbreg:::backend_get_query(remote, paste0("SELECT ", dbreg:::sql_stable_sum(remote, "v"),
      " AS score,", dbreg:::sql_sum(remote, "v"), " AS expected FROM (VALUES ", values, ") t(v)"))
    expect_true(is.finite(row[["score"]]) && row[["score"]] != 0 &&
      abs((row[["score"]] - row[["expected"]]) / row[["expected"]]) < 1e-12)
  }
})
