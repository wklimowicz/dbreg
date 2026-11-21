# A rare high-risk group can make the first Newton step improve total
# deviance while clamping its two observed failures to probability one.
local({
  set.seed(271)
  n = 4000L
  d = data.frame(g = rep(1:40, length.out = n), h = sample(1:8, n, TRUE), x = rnorm(n))
  d[["y"]] = rbinom(n, 1, -expm1(-exp(-4.5 + 0.2 * d[["x"]] + 0.1 * sin(d[["h"]]))))
  d[["y"]][d[["g"]] == 1L] = c(0L, 0L, rep(1L, 98L))
  family = binomial("cloglog")
  # Match perfect-fit pruning in the reference by dropping zero-event groups.
  keep = d[["g"]] %in% d[["g"]][d[["y"]] == 1L]
  reference = stats::glm(y ~ x + factor(g) + factor(h), data = d[keep, ], family = family,
    control = stats::glm.control(epsilon = 1e-14))
  expect_true(reference[["converged"]])
  conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))

  for (method in c("auto", "duckdb_bounded")) {
    fit = dbglm(y ~ x | g + h, conn = conn, data = d, family = family,
      control = dbglm_control(projection_method = method, score_tol = 1e-9))
    expect_true(fit[["converged"]] && fit[["finite_mle"]])
    expect_true(fit[["step_halvings"]] > 0L)
    expect_true(max(fit[["score_max"]], fit[["fe_score_max"]]) <= 1e-9)
    expect_equal(nobs(fit), nobs(reference))
    expect_equal(coef(fit), coef(reference)["x"], tolerance = 1e-6)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)), tolerance = 1e-8)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})
