# A few tiny row weights need not make a well-supported additive design weak.
local({
  set.seed(1751)
  n = 2400L
  d = data.frame(
    g = sample(16L, n, TRUE), h = sample(9L, n, TRUE),
    district = sample(12L, n, TRUE), x = rnorm(n), w = 1
  )
  d[["w"]][seq.int(1L, n, 40L)] = 1e-12
  eta = -1.4 + .3 * d[["x"]] + .15 * sin(d[["g"]]) + .1 * cos(d[["h"]])
  d[["y"]] = rbinom(n, 1L, -expm1(-exp(eta)))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))

  fit = function(method) {
    dbglm(y ~ x | g + h, conn = conn, data = d, weights = ~w,
      family = binomial("cloglog"), vcov = ~district, strategy = "moments",
      control = dbglm_control(projection_method = method, score_tol = 1e-9))
  }
  reference = fit("auto")
  candidate = fit("duckdb_bounded")
  expect_identical(reference[["projection_diagnostics"]][["reason"]], "weak_information")
  expect_identical(candidate[["projection_method"]], "bounded_moments")
  expect_true(candidate[["converged"]] && candidate[["finite_mle"]])
  expect_true(max(candidate[["score_max"]], candidate[["fe_score_max"]]) <= 1e-9)
  expect_equal(nobs(candidate), nobs(reference))
  expect_equal(coef(candidate), coef(reference), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(candidate)), as.numeric(vcov(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(candidate)), as.numeric(logLik(reference)), tolerance = 1e-8)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
