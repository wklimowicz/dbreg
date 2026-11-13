# dbglm accepts a lazy table as `data`, as dbreg does.
if (!requireNamespace("dbplyr", quietly = TRUE)) exit_file("dbplyr not installed")
local({
  set.seed(3); n = 800
  d = data.frame(id = 1:n, x = rnorm(n), g = sample(20, n, TRUE))
  d$y = rpois(n, exp(.3 * d$x + d$g / 40))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  lazy = dplyr::filter(dplyr::tbl(conn, "d"), x > -1)
  ref = dbglm(y ~ x | g, conn = conn, table = lazy, family = poisson())
  fit = dbglm(y ~ x | g, data = lazy, family = poisson())
  expect_equal(coef(fit), coef(ref))
  expect_equal(nobs(fit), sum(d$x > -1))
  fit = dbglm(y ~ x | g, data = lazy, family = poisson(), row_id = "id", sample_table = "glm_lazy_sample")
  expect_equal(nrow(DBI::dbReadTable(conn, "glm_lazy_sample")), sum(d$x > -1))
})
