# Weakly connected teacher/school effects converge slowly with plain AP.
# Acceleration must meet the same certificate within the same cycle budget.
local({
  saved = options(dbreg.ap_max_iter = 200L, dbreg.ap_tol = 1e-10,
    dbreg.ap_method = getOption("dbreg.ap_method"),
    dbreg.ap_accelerate = getOption("dbreg.ap_accelerate"))
  on.exit(options(saved), add = TRUE)
  set.seed(1)
  nt = 300L; years = 6L
  d = data.frame(tid = rep(seq_len(nt), each = years), year = rep(seq_len(years), nt))
  home = sample.int(400, nt, TRUE)
  move = runif(nrow(d)) < .03
  d[["school"]] = ifelse(move, sample.int(400, nrow(d), TRUE), home[d[["tid"]]])
  d[["age"]] = pmin(60, 25 + (d[["tid"]] %% 30) + d[["year"]])
  d[["x"]] = rnorm(nrow(d)) + .3*(d[["school"]] %% 7) + .1*d[["age"]]
  d[["y"]] = .5*d[["x"]] + rnorm(nt)[d[["tid"]]] + rnorm(400)[d[["school"]]] + rnorm(nrow(d))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "teacher_school", d)
  fml = y ~ x | tid + school + year + age
  ref = fixest::feols(fml, d, fixef.tol = 1e-11, fixef.iter = 100000L, fixef.rm = "none")
  for (method in c("join", "window")) {
    options(dbreg.ap_method = method, dbreg.ap_accelerate = FALSE)
    expect_error(dbreg(fml, conn = conn, table = "teacher_school", strategy = "demean"), "did not converge")
    expect_identical(DBI::dbListTables(conn), "teacher_school")
    options(dbreg.ap_accelerate = TRUE)
    fit = dbreg(fml, conn = conn, table = "teacher_school", strategy = "demean")
    expect_true(fit[["ap_iterations"]] <= 200L)
    expect_true(fit[["ap_accelerated_steps"]] > 0L)
    expect_true(fit[["ap_max_abs"]] < 1e-10)
    expect_true(fit[["ap_max_relative"]] <= 1e-10)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-9)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
    expect_identical(DBI::dbListTables(conn), "teacher_school")
  }
  for (bad in list("yes", NA, logical(), c(TRUE, FALSE))) {
    options(dbreg.ap_accelerate = bad)
    expect_error(dbreg(fml, conn = conn, table = "teacher_school", strategy = "demean"), "dbreg.ap_accelerate")
    expect_identical(DBI::dbListTables(conn), "teacher_school")
  }
})
