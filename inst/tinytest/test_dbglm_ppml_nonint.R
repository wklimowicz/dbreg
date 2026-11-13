# PPML accepts nonnegative non-integer outcomes (trade values), as fepois and glm do.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(8)
  d = expand.grid(ex = 1:25, im = 1:25)
  d = d[d$ex != d$im, ]
  d$ldist = log(runif(nrow(d), 100, 1e4))
  mu = exp(2 + rnorm(25)[d$ex] + rnorm(25)[d$im] - 0.75 * d$ldist + 8)
  d$flow = round(rpois(nrow(d), mu) * runif(nrow(d), 0.9, 1.1), 3)
  d$flow[sample(nrow(d), 30)] = 0
  ref = fixest::fepois(flow ~ ldist | ex + im, d, glm.tol = 1e-12, fixef.tol = 1e-11, notes = FALSE)
  fit = dbglm(flow ~ ldist | ex + im, data = d, family = poisson())
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(sqrt(diag(vcov(fit)))), as.numeric(fixest::se(ref)), tolerance = 1e-6)
  expect_equal(deviance(fit), deviance(ref), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-8)
  # The portable (no lgamma builtin) saturated term agrees with lgamma.
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  y = c(0, 0.3, 1, 2.5, 7, 19.99, 20, 55.5, 1000.25)
  DBI::dbWriteTable(conn, "b", data.frame(`__dbglm_y` = y, `__dbglm_a` = 1, check.names = FALSE))
  portable = DBI::dbGetQuery(conn, paste("SELECT", dbreg:::dbglm_saturated_loglik_term(FALSE), "AS t FROM b"))[["t"]]
  expect_equal(portable, ifelse(y == 0, 0, y * log(y) - y - lgamma(y + 1)), tolerance = 1e-12)
})
