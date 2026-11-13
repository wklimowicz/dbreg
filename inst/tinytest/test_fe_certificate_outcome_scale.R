# Fixed-effect certificates are relative to the column scale, so dollar-scale
# outcomes certify on the bounded and coordinate engines.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(21); n = 3000
  d = data.frame(x0 = rnorm(n), z = rnorm(n), t = rnorm(n), firm = sample(200, n, TRUE, prob = rexp(200)), year = sample(30, n, TRUE))
  e0 = 0.5 * d$x0 - 0.2 * d$z + rnorm(200)[d$firm] + rnorm(30)[d$year] + 0.4 * d$t * rnorm(200)[d$firm] + rnorm(n)
  d$y = 1e8 * (1 + 0.1 * e0)
  d$x = 1e6 * (5 + d$x0)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  for (f in list(y ~ x + z | firm + year, y ~ x + z | firm[t] + year)) {
    ref = fixest::feols(f, d, fixef.tol = 1e-11, fixef.iter = 1e5)
    fit = dbreg(f, conn = conn, table = "d")
    expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8, info = deparse(f))
  }
  fit = dbreg(y ~ x + z | firm + year, conn = conn, table = "d")
  expect_false(grepl("declined", paste(unlist(fit[["projection_diagnostics"]]), collapse = " ")))
})
