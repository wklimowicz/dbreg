source(system.file("tinytest/helper_ap_moments.R", package = "dbreg"), local = TRUE)
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_ap_moments(conn)
  expect_equal(length(DBI::dbListTables(conn)), 0L)

  # Exercise ownership when the shared solver cannot reconstruct its result.
  # A cloned function environment injects failures without changing a namespace.
  DBI::dbWriteTable(conn, "source", data.frame(f = rep(1:3, each = 4),
    g = rep(1:2, 6), x = seq_len(12), y = sin(seq_len(12))))
  inputs = list2env(list(conn = conn, from_statement = "FROM source", fe = c("f", "g"),
    yvar = "y", weights = NULL, cluster_var = NULL, verbose = FALSE))
  project = dbreg:::dbreg_moment_projections
  environment(project) = new.env(parent = environment(project))
  injected = "FE coefficient lookup exceeds the SQL transfer limit."
  environment(project)[["project_fe_moments"]] = function(...) dbreg:::projection_abort("upload_limit", injected)
  attempt = dbreg:::projection_attempt("bounded_moments", project(inputs, "x", "x", 100L, 1e-10))
  expect_null(attempt[["result"]])
  expect_identical(attempt[["reason"]], "upload_limit")
  expect_identical(DBI::dbListTables(conn), "source")
  injected = "injected unexpected database failure"
  environment(project)[["project_fe_moments"]] = function(...) stop(injected)
  expect_error(project(inputs, "x", "x", 100L, 1e-10), injected)
  expect_identical(DBI::dbListTables(conn), "source")
  environment(project)[["detect_backend"]] = function(conn) list(name = "spark",
    is_databricks = TRUE, supports_physical_temp_tables = TRUE, native_query = FALSE)
  expect_error(project(inputs, "x", "x", 100L, 1e-10), "native query mode")
  expect_identical(DBI::dbListTables(conn), "source")
})
