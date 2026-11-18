# Each SQL CG cycle materializes a single state table (replaced in place on
# DuckDB); certificates and scalar diagnostics never create tables.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  counts = c(create = 0L, replace = 0L, drop = 0L)
  create = function(sql) { counts[["create"]] <<- counts[["create"]] + 1L; manager$create(sql) }
  drop = function(table) { counts[["drop"]] <<- counts[["drop"]] + 1L; manager$drop(table) }
  replace = function(table, sql) {
    counts[[if (is.null(table)) "create" else "replace"]] <<- counts[[if (is.null(table)) "create" else "replace"]] + 1L
    manager$replace(table, sql)
  }
  set.seed(11147)
  d = expand.grid(g = 1:6, h = 1:5, k = 1:4)
  d[["w"]] = runif(nrow(d), .25, 2)
  d[["x"]] = rnorm(nrow(d)) + d[["g"]]
  d[["row_id"]] = seq_len(nrow(d))
  DBI::dbWriteTable(conn, "cg_count", d)
  fe = c("g", "h", "k")
  X = model.matrix(~ factor(g) + factor(h) + factor(k), d)
  oracle = lm.wfit(X, d[["x"]], d[["w"]])[["residuals"]]
  fit = dbreg:::fe_sql_cg(conn, "cg_count", "x", c("row_id", fe),
    lapply(fe, function(f) list(f = f, v = "1.0")), 1, 1e-12,
    dbreg:::fe_cg_control(), create, drop,
    list2env(list(rows = 0, cells = 0), parent = emptyenv()),
    dbreg:::projection_limits(), retain_coefficients = TRUE, replace = replace)
  expect_true(fit[["cycles"]] > 2L)
  # Group scales, the first state, the second ping-pong state, projected rows,
  # and coefficients are created once; every later cycle replaces a state.
  # A restart (also a stalled final step) adds one certified residual table.
  extra = counts[["create"]] - 5L
  expect_true(extra >= 0L && extra <= fit[["restarts"]] + 1L)
  expect_equal(counts[["replace"]], fit[["cycles"]] - 2L)
  expect_equal(counts[["drop"]], 3L + extra)
  residuals = DBI::dbGetQuery(conn, paste("SELECT x FROM", fit[["table"]], "ORDER BY row_id"))[["x"]]
  expect_equal(residuals, unname(oracle), tolerance = 1e-10)
  # Warm-starting from the converged coefficients needs (almost) no cycles.
  warm = dbreg:::fe_sql_cg(conn, "cg_count", "x", c("row_id", fe),
    lapply(fe, function(f) list(f = f, v = "1.0")), 1, 1e-12,
    dbreg:::fe_cg_control(), create, drop,
    list2env(list(rows = 0, cells = 0), parent = emptyenv()),
    dbreg:::projection_limits(), replace = replace,
    start = list(table = fit[["coefficients"]], ratio = 1))
  expect_true(warm[["cycles"]] < fit[["cycles"]])
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT x FROM", warm[["table"]], "ORDER BY row_id"))[["x"]],
    unname(oracle), tolerance = 1e-10)
  for (table in c(fit[["coefficients"]], fit[["table"]], warm[["table"]])) manager$drop(table)
  expect_identical(DBI::dbListTables(conn), "cg_count")
})
