# Explicit compressed varying-slope fits certify residual drift relative to the outcome scale.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(12); n = 6000
  d = data.frame(id = sample(50, n, TRUE), firm = sample(20, n, TRUE), age = sample(20:60, n, TRUE),
    hours = sample(c(20, 30, 40, 50), n, TRUE))
  d$wage = 15000 + 50 * d$age + rnorm(50, 0, 5000)[d$id] + rnorm(20, 0, 2000)[d$firm] +
    rnorm(20, 300, 50)[d$firm] * d$hours + rnorm(n, 0, 3000)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  ref = fixest::feols(wage ~ age | id + firm[hours], d, fixef.tol = 1e-11, notes = FALSE)
  fit = dbreg(wage ~ age | id + firm[hours], conn = conn, table = "d", strategy = "compress")
  expect_equal(fit[["strategy"]], "compress")
  expect_equal(unname(coef(fit)["age"]), unname(coef(ref)["age"]), tolerance = 1e-8)
})

# DECIMAL money regressors, revenue at 1e8 and clustered SEs on compressed varying slopes.
local({
  set.seed(11); n = 4000
  d = data.frame(firm = sample(200, n, TRUE), yr = sample(2010:2019, n, TRUE), age = sample(18:65, n, TRUE),
    cl = sample(letters, n, TRUE))
  d$price = round(runif(n, 1, 500), 2)
  d$rev = round(1e8 + 3000 * d$price + 50 * d$age * rnorm(200, 1, .2)[d$firm] + rnorm(200, 0, 2e5)[d$firm] +
    rnorm(10, 0, 1e5)[d$yr - 2009] + rnorm(n, 0, 3e5), 2)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "draw", d)
  DBI::dbExecute(conn, "CREATE TABLE d AS SELECT firm, yr, age, cl, CAST(price AS DECIMAL(10,2)) AS price, rev FROM draw")
  for (V in list("iid", ~cl)) {
    ref = fixest::feols(rev ~ price | firm[age] + yr, d, vcov = V, fixef.tol = 1e-11, ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
    fit = dbreg(rev ~ price | firm[age] + yr, conn = conn, table = "d", strategy = "compress", vcov = V)
    expect_equal(fit[["strategy"]], "compress")
    expect_equal(unname(coef(fit)["price"]), unname(coef(ref)["price"]), tolerance = 1e-8)
    expect_equal(unname(sqrt(diag(vcov(fit)))[1]), unname(fixest::se(ref)["price"]), tolerance = 1e-7)
  }
})
