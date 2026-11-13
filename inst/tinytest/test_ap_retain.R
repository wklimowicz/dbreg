# Retained OLS fits keep their row keys through alternating projections.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(4); n = 2000
  d = data.frame(id = seq_len(n), x = rnorm(n), worker = sample(300, n, TRUE), year = sample(2001:2020, n, TRUE))
  d$y = .5 * d$x + rnorm(300)[d$worker] + rnorm(20)[d$year - 2000] + rnorm(n)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  ref = fixest::feols(y ~ x | worker + year, d, fixef.tol = 1e-11, ssc = fixest::ssc(K.fixef = "full"))
  old = options(dbreg.ap_solver = "iterative"); on.exit(options(old), add = TRUE)
  for (retain in c("prediction", "inference", "all")) {
    fit = dbreg(y ~ x | worker + year, conn = conn, table = "d", retain = retain, row_id = "id")
    expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8, info = retain)
    if (retain != "inference") expect_equal(as.numeric(predict(fit, newdata = d[1:5, ])), as.numeric(predict(ref, newdata = d[1:5, ])),
      tolerance = 1e-6, info = retain)
    close(fit)
  }
})
