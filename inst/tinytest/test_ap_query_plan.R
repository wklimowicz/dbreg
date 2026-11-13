# Plans create no AP tables. Their ordinary-AP templates execute to the same
# certified fit, including the original-scale rank guard for an FE-only term.
local({
  saved = options(dbreg.ap_method = getOption("dbreg.ap_method"))
  on.exit(options(saved), add = TRUE)
  set.seed(850)
  d = data.frame(x = rnorm(150), z = rnorm(150), f = sample(1:7, 150, TRUE),
    g = sample(1:5, 150, TRUE), w = runif(150, .4, 2))
  d[["absorbed"]] = sin(d[["f"]]) + cos(d[["g"]])
  d[["y"]] = d[["x"]] + d[["absorbed"]] + rnorm(150)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "plan_fixture", d)
  ref = fixest::feols(y ~ x + z | f + g, d, weights = ~w, fixef.rm = "none")
  for (method in c("join", "window")) {
    options(dbreg.ap_method = method)
    plan = dbreg(y ~ x + z + absorbed | f + g, conn = conn, table = "plan_fixture",
      weights = "w", vcov = ~x, strategy = "demean", sql_only = TRUE, fe_tol = 1e-9)
    expect_true(inherits(plan, "dbreg_query_plan"))
    expect_identical(DBI::dbListTables(conn), "plan_fixture")
    expect_identical(plan[["method"]], method)
    expect_true(length(capture.output(print(plan))) > 0L)
    init = gsub("<model_frame_sql>", plan[["model_frame_sql"]], plan[["initialization_sql"]], fixed = TRUE)
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE plan_current AS", init))
    scales = as.numeric(DBI::dbGetQuery(conn, gsub("<model_frame_sql>", plan[["model_frame_sql"]], plan[["original_scale_sql"]], fixed = TRUE))[1, ])
    for (iter in seq_len(plan[["max_iterations"]])) {
      for (projection in plan[["projection_sql"]]) {
        sql = gsub("<ap_current>", "plan_current", projection, fixed = TRUE)
        DBI::dbExecute(conn, paste("CREATE OR REPLACE TEMP TABLE plan_next AS", sql))
        DBI::dbExecute(conn, "CREATE OR REPLACE TEMP TABLE plan_current AS SELECT * FROM plan_next")
      }
      absolute = as.numeric(DBI::dbGetQuery(conn, gsub("<ap_current>", "plan_current", plan[["convergence_sql"]], fixed = TRUE))[1, ])
      rms = as.numeric(DBI::dbGetQuery(conn, gsub("<ap_current>", "plan_current", plan[["residual_scale_sql"]], fixed = TRUE))[1, ])
      identified = rms > 128*.Machine[["double.eps"]]*scales
      relative = if (any(identified)) max(absolute[identified]/rms[identified]) else 0
      if (max(absolute) < plan[["tolerance"]] && relative <= plan[["tolerance"]]) break
    }
    expect_true(iter < plan[["max_iterations"]])
    moments = DBI::dbGetQuery(conn, gsub("<ap_converged>", "plan_current", plan[["final_moments_sql"]], fixed = TRUE))
    S = dbreg:::cross_moments(moments, c("x", "z", "absorbed"), "y", intercept = FALSE)
    fit = dbreg:::fit_ols(S, "y", c("x", "z", "absorbed"), drop_collinear = TRUE)
    expect_identical(fit[["dropped"]], "absorbed")
    expect_equal(as.numeric(fit[["coef"]]), as.numeric(coef(ref)), tolerance = 1e-8)
    DBI::dbExecute(conn, "DROP TABLE plan_current")
    DBI::dbExecute(conn, "DROP TABLE plan_next")
  }
  iv = dbivreg(y ~ x | f + g | z ~ w, conn = conn, table = "plan_fixture", sql_only = TRUE)
  expect_true(inherits(iv, "dbreg_query_plan"))
  expect_true(all(c("x", "z", "w") %in% iv[["design_matrix"]][["alias"]]))
  expect_identical(DBI::dbListTables(conn), "plan_fixture")
  ordinary = capture.output({
    sql = dbreg(y ~ x, conn = conn, table = "plan_fixture", strategy = "moments", sql_only = TRUE)
  })
  expect_true(is.character(sql))
  expect_true(length(ordinary) > 0L)
})
