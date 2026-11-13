# augment() returns the fitted values and residuals of the same fixest model.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(20261117)
  n = 600
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = sample(25, n, TRUE),
    h = sample(8, n, TRUE), w = runif(n, .5, 2))
  d$y = .5*d$x - .2*d$z + rnorm(25)[d$g] + rnorm(8)[d$h] + rnorm(n)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  for (f in list(y ~ x + z | g, y ~ x + z | g + h)) for (s in c("auto", "compress")) {
    ref = fixest::feols(f, d, weights = ~w, fixef.tol = 1e-11, fixef.rm = "none")
    fit = dbreg(f, conn = conn, table = "d", weights = "w", strategy = s, retain = "inference", row_id = "id")
    a = generics::augment(fit); a = a[order(a$id), ]
    expect_equal(a$.fitted, unname(fitted(ref)), tolerance = 1e-8, info = paste(deparse(f), s))
    expect_equal(a$.resid, unname(resid(ref)), tolerance = 1e-7, info = paste(deparse(f), s))
    close(fit)
  }
})

## Compressed fits that retain inference rows fall back to demean when the cell scores
## would exceed the SQL upload budget.
local({
  set.seed(1)
  n = 4000
  d = data.frame(rid = 1:n, x = rnorm(n), z = rnorm(n), g = sample(200, n, TRUE))
  d$y = .3 * d$x - .2 * d$z + rnorm(200)[d$g] + rnorm(n)
  con = DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "d", d)
  before = DBI::dbListTables(con)
  ref = fixest::feols(y ~ x + z | g, d, vcov = ~g, ssc = fixest::ssc(K.fixef = "full"))
  fit = dbreg(y ~ x + z | g, conn = con, table = "d", strategy = "compress", vcov = ~g,
    retain = "inference", row_id = "rid")
  expect_equal(fit[["strategy"]], "demean")
  expect_equal(coef(fit)[c("x", "z")], coef(ref)[c("x", "z")], tolerance = 1e-10)
  expect_equal(sqrt(diag(vcov(fit)))[c("x", "z")], fixest::se(ref)[c("x", "z")], tolerance = 1e-10)
  res = residuals(fit)
  expect_equal(as.data.frame(res)[[2L]], as.numeric(residuals(ref)), tolerance = 1e-8)
  close(res)
  close(fit)
  expect_identical(DBI::dbListTables(con), before)
  DBI::dbDisconnect(con, shutdown = TRUE)
})
