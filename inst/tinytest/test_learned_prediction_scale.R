# Learned-effect prediction works for dollar-scale outcomes, factor-by-factor
# cells beside a FE, and three fixed-effect families.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(3); n = 4000
  d = data.frame(id = sample(200, n, TRUE), edu = sample(c("hs", "ba", "ma"), n, TRUE),
    region = sample(c("N", "S", "E"), n, TRUE), yr = sample(8, n, TRUE), hours = sample(20:60, n, TRUE))
  d$wage = 40000 + 900 * d$hours + rnorm(200, 0, 15000)[d$id] + 2000 * d$yr + rnorm(n, 0, 10000)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  nd = d[1:20, ]
  for (f in list(wage ~ hours + edu:region | id, wage ~ hours | id + yr + region)) {
    ref = fixest::feols(f, d, fixef.tol = 1e-11, notes = FALSE)
    expect_silent(fit <- dbreg(f, conn = conn, table = "d", retain = "prediction"))
    expect_equal(as.numeric(predict(fit, newdata = nd)), as.numeric(predict(ref, newdata = nd)),
      tolerance = 1e-8, info = deparse(f))
    close(fit)
  }
})

# Learned-effect dictionaries for worker trends and three FE families fit the
# default prediction budget (it counts dictionary entries, not a dense square).
local({
  set.seed(14)
  d = expand.grid(worker = 1:1200, year = 1:5)
  d$tenure = d$year + sample(0:10, 1200, TRUE)[d$worker]
  d$firm = sample(300, nrow(d), TRUE)
  d$hours = rnorm(nrow(d), 40, 5)
  d$wage = 20 + 0.1 * d$hours + rnorm(1200)[d$worker] + rnorm(1200, 0, .1)[d$worker] * d$tenure +
    rnorm(300)[d$firm] + 0.3 * d$year + rnorm(nrow(d))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  nd = d[1:20, ]
  for (f in list(wage ~ hours | worker[tenure] + year, wage ~ hours | worker + firm + year)) {
    ref = fixest::feols(f, d, fixef.tol = 1e-11, notes = FALSE)
    expect_silent(fit <- dbreg(f, conn = conn, table = "d", retain = "prediction"))
    expect_equal(as.numeric(predict(fit, newdata = nd)), as.numeric(predict(ref, newdata = nd)),
      tolerance = 1e-7, info = deparse(f))
    close(fit)
  }
})
