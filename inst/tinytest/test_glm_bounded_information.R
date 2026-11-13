# Bound-based guard elision must preserve actual guarded aggregate values.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  native = list(name = "spark", is_databricks = TRUE)
  sum_sql = function(bound = NULL) dbreg:::sql_stable_sum(conn, "v", native, max_abs = bound)
  expected = c(1e-4, 2e-27, 625)
  index = 0L
  for (values in list(c(1e5, 1e-7, -1e5), c(1e-30, -1e-30, 2e-30), c(-.125, .5, .25))) {
    DBI::dbWriteTable(conn, "numbers", data.frame(v = rep(values, 1000)))
    index = index + 1L
    bound = max(abs(values))
    ordinary = DBI::dbGetQuery(conn, paste("SELECT", sum_sql(), "AS total FROM numbers"))$total
    fast = DBI::dbGetQuery(conn, paste("SELECT", sum_sql(bound), "AS total FROM numbers"))$total
    expect_equal(fast, ordinary, tolerance = 1e-14)
    expect_equal(fast, expected[index], tolerance = 1e-12)
    expect_true(nchar(sum_sql(bound)) < nchar(sum_sql()))
    DBI::dbRemoveTable(conn, "numbers")
  }
  DBI::dbWriteTable(conn, "moments", data.frame(w = 1,
    a = rep(c(1e5, 1e-7, -1e5), 1000), b = 1, tiny = 1e-30))
  query = dbreg:::fe_information_query
  qe = new.env(parent = environment(query)); qe$detect_backend = function(...) native
  environment(query) = qe
  fresh = function() list2env(list(rows = 0, cells = 0))
  unbounded = query(conn, "moments", c("a", "b", "tiny"), fresh(), dbreg:::projection_limits())
  bounded = query(conn, "moments", c("a", "b", "tiny"), fresh(), dbreg:::projection_limits(),
    list(weight = 1, columns = c(1e5, 1, 1e-30)))
  expect_equal(bounded, unbounded, tolerance = 1e-12)
  expect_equal(bounded$cross[1, 2], 1e-4, tolerance = 1e-12)
  expect_true(bounded$cross[3, 3] > 0)
  expect_equal(bounded$cross[3, 3] / 1e-60, 3000, tolerance = 1e-12)
  DBI::dbRemoveTable(conn, "moments")
  original = sum_sql()
  for (bound in list(NULL, NA_real_, Inf, NaN, -1, 1e9, 1e10, "invalid")) expect_identical(sum_sql(bound), original)
  expect_true(nchar(sum_sql(1e9 - 1)) < nchar(original))
  expect_identical(dbreg:::sql_stable_sum(conn, "v", max_abs = 1), dbreg:::sql_stable_sum(conn, "v"))
})

local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = expand.grid(g = 1:2, h = 1:2, x = c(-1, 1), rep = 1:2)
  d$w = 1; d$score = sin(seq_len(nrow(d)))
  DBI::dbWriteTable(conn, "working", d)
  step = dbreg:::dbglm_fe_step
  e = new.env(parent = environment(step))
  e$detect_backend = function(...) list(name = "spark", is_databricks = TRUE)
  capture = new.env(parent = emptyenv())
  e$fe_fwl_information = function(conn, source, columns, B, coefficients, sys, budget, limits, global = NULL, bounds = NULL) {
    capture$bounds = bounds
    dbreg:::fe_fwl_information(conn, source, columns, B, coefficients, sys, budget, limits, global, bounds)
  }
  environment(step) = e
  fit = step(conn, "working", "x", c("g", "h"), c(2L, 2L), nrow(d), 1e-9)
  expect_equal(capture$bounds$weight, 1)
  expect_equal(capture$bounds$columns, c(max(abs(d$score)), 1))
  expect_true(fit$certificate$max_abs < 1e-9)
  for (field in c("w", "score", "x")) for (bad in c(NA_real_, NaN, Inf)) {
    invalid = d; invalid[[field]][1L] = bad
    DBI::dbWriteTable(conn, "invalid", invalid, overwrite = TRUE)
    error = tryCatch(step(conn, "invalid", "x", c("g", "h"), c(2L, 2L), nrow(d), 1e-9), error = identity)
    expect_true(inherits(error, "dbreg_projection_numerical"))
    DBI::dbRemoveTable(conn, "invalid")
  }
  limited = step
  le = new.env(parent = e)
  le$projection_limits = function() {
    limits = dbreg:::projection_limits(); limits$transfer_cells = 4; limits
  }
  environment(limited) = le
  error = tryCatch(limited(conn, "must_not_query", "x", c("g", "h"),
    c(2L, 2L), nrow(d), 1e-9), error = identity)
  expect_true(inherits(error, "dbreg_projection_error"))
  expect_identical(error$reason, "transfer_limit")
  # Cached Fisher states retain the original four-cell norm read and reach SQL.
  cached = list(working = "must_not_query", working_columns = "x",
    information_type = "fisher", H = matrix(1, 1, 1), score = 0)
  error = tryCatch(limited(conn, "must_not_query", "x", c("g", "h"),
    c(2L, 2L), nrow(d), 1e-9, raw_state = cached), error = identity)
  expect_true(inherits(error, "error"))
  expect_false(inherits(error, "dbreg_projection_numerical"))
  # Projected coordinates have different bounds; actual-H fallback is unbounded.
  info = dbreg:::dbglm_fe_information
  ie = new.env(parent = environment(info))
  ie$fe_information_query = function(conn, source, columns, budget, limits, bounds = NULL) {
    capture$fallback_bounds = bounds
    dbreg:::fe_information_query(conn, source, columns, budget, limits, bounds)
  }
  environment(info) = ie
  DBI::dbExecute(conn, "CREATE TEMP TABLE residual AS SELECT w,score AS v1,x AS v2 FROM working")
  value = info(conn, "residual", c("v1", "v2"),
    list(cross = diag(2), within = diag(c(1, -1))),
    list2env(list(cells = 0, rows = 0)), dbreg:::projection_limits(), nrow(d))
  expect_true(value$recomputed)
  expect_null(capture$fallback_bounds)
  DBI::dbRemoveTable(conn, "residual")
  DBI::dbRemoveTable(conn, "working")
  expect_identical(DBI::dbListTables(conn), character())
})
