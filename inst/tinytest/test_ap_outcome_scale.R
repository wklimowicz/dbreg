# Alternating projections certify group means relative to the column scale, so
# dollar-scale outcomes converge on every backend that uses them.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(5); n = 3000
  d = data.frame(x = rnorm(n), firm = sample(150, n, TRUE, prob = rexp(150)), yr = sample(20, n, TRUE),
    state = sample(12, n, TRUE), w = runif(n, 1, 1000))
  d$y = 3e8 * (1 + 0.05 * (0.5 * d$x + rnorm(150)[d$firm] + rnorm(20)[d$yr] + rnorm(n)))
  ref = fixest::feols(y ~ x | firm + yr, d, fixef.tol = 1e-11, fixef.iter = 1e5)
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    conn = DBI::dbConnect(RSQLite::SQLite(), ":memory:")
    DBI::dbWriteTable(conn, "d", d)
    fit = dbreg(y ~ x | firm + yr, conn = conn, table = "d", strategy = "demean")
    expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
    DBI::dbDisconnect(conn)
  }
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  old = options(dbreg.ap_solver = "iterative"); on.exit(options(old), add = TRUE)
  ref3 = fixest::feols(y ~ x | firm + yr + state, d, weights = ~w, fixef.tol = 1e-11, fixef.iter = 1e5)
  fit3 = dbreg(y ~ x | firm + yr + state, conn = conn, table = "d", weights = "w", strategy = "demean")
  expect_equal(unname(coef(fit3)), unname(coef(ref3)), tolerance = 1e-8)
})
