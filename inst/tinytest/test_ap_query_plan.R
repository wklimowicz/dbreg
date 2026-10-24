set.seed(850)
dat = data.frame(x = rnorm(150), z = rnorm(150), f = sample(1:7, 150, TRUE),
  g = sample(1:5, 150, TRUE), w = runif(150, 0.4, 2))
dat[["y"]] = dat[["x"]] + dat[["f"]] + rnorm(150)
con = DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con, "fixture", dat)
old = options(dbreg.ap_solver = "iterative", dbreg.ap_method = getOption("dbreg.ap_method"))
for (method in c("join", "window")) {
  options(dbreg.ap_method = method)
  plan = dbreg(y ~ x + z | f + g, conn = con, table = "fixture", weights = "w",
    cluster = ~x, strategy = "demean", sql_only = TRUE, fe_tol = 1e-9)
  expect_true(inherits(plan, "dbreg_query_plan"))
  expect_equal(DBI::dbListTables(con), "fixture")
  init = gsub("<model_frame_sql>", plan[["model_frame_sql"]], plan[["initialization_sql"]], fixed = TRUE)
  DBI::dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE review_ap AS", init))
  for (iter in seq_len(100L)) {
    for (projection in plan[["projection_sql"]]) {
      sql = gsub("<ap_current>", "review_ap", projection, fixed = TRUE)
      DBI::dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE review_next AS", sql))
      DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE review_ap AS SELECT * FROM review_next")
    }
    convergence = DBI::dbGetQuery(con, gsub("<ap_current>", "review_ap", plan[["convergence_sql"]], fixed = TRUE))
    if (max(as.numeric(convergence[1, ])) < plan[["tolerance"]]) break
  }
  expect_true(iter < 100L)
  moments = DBI::dbGetQuery(con, gsub("<ap_converged>", "review_ap", plan[["final_moments_sql"]], fixed = TRUE))
  S = dbreg:::cross_moments(moments, c("x", "z"), "y", intercept = FALSE)
  beta = drop(solve(S[c("x", "z"), c("x", "z")], S[c("x", "z"), "y"]))
  ref = fixest::feols(y ~ x + z | f + g, dat, weights = ~w, fixef.rm = "none")
  expect_equal(beta, coef(ref), tolerance = 1e-8)
  DBI::dbExecute(con, "DROP TABLE review_ap")
  DBI::dbExecute(con, "DROP TABLE review_next")
}
plan_iv = dbivreg(y ~ x | f + g | z ~ w, conn = con, table = "fixture", sql_only = TRUE)
expect_true(inherits(plan_iv, "dbreg_query_plan"))
expect_equal(DBI::dbListTables(con), "fixture")
options(old)

DBI::dbDisconnect(con, shutdown = TRUE)
