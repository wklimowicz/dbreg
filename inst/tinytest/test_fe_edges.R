source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_fe_edges.R", package = "dbreg"), local = TRUE)
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_fe_edges(conn, "glm")
  dbreg_test_fe_edges(conn, "ols")
  # nobs() counts source observations, never moment rows or compressed cells.
  d = expand.grid(x = 1:5, f = 1:4, repeat_id = 1:3)
  d[["y"]] = .4 * d[["x"]] + sin(d[["f"]]) + cos(d[["repeat_id"]])
  for (strategy in c("moments", "compress", "demean", "mundlak")) {
    formula = if (strategy == "moments") y ~ x else y ~ x | f
    got = dbreg(formula, conn = conn, data = d, strategy = strategy)
    expect_equal(nobs(got), nrow(d), info = strategy)
    expect_equal(nobs(got), generics::glance(got)[["nobs"]], info = strategy)
  }
  expect_equal(length(DBI::dbListTables(conn)), 0L)
})
