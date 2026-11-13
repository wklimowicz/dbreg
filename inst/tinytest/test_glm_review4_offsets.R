# Optional null diagnostics must reuse the offset-aware initializer that lets
# the main fit recover from a saturated mean-offset default.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(y = c(rep(0, 20), rep(c(0, 1), 10)),
    o = rep(c(-100, 0), each = 20), g = rep(rep(1:2, each = 5), 4))
  for (family in list(binomial(), binomial("probit"), binomial("cloglog"))) {
    for (strategy in c("moments", "compress")) {
      fit = dbglm(y ~ 1 | g, conn = conn, data = d, offset = "o", family = family,
        strategy = strategy, null_model = TRUE)
      null = dbglm(y ~ 1, conn = conn, data = d, offset = "o", family = family,
        strategy = strategy)
      expect_true(fit[["null_diagnostics"]][["fit"]][["available"]])
      expect_equal(fit[["null_diagnostics"]][["weighted"]][["logLik"]],
        as.numeric(logLik(null, type = "weighted")), tolerance = 1e-8)
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
  }
})

# The Poisson constant score has a stable closed-form solution with offsets.
# A mean-offset start instead needs about 300 unit scoring steps here.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(y = rep(c(1, 2), 20), o = rep(c(-600, 0), each = 20),
    g = rep(rep(1:2, each = 5), 4))
  expected = log(sum(d[["y"]]) / sum(exp(d[["o"]])))
  for (strategy in c("moments", "compress")) {
    fit = dbglm(y ~ 1, conn = conn, data = d, offset = "o", family = poisson(), strategy = strategy)
    expect_true(fit[["converged"]] && fit[["finite_mle"]])
    expect_equal(unname(coef(fit)), expected, tolerance = 1e-10)
    absorbed = dbglm(y ~ 1 | g, conn = conn, data = d, offset = "o", family = poisson(),
      strategy = strategy, null_model = TRUE)
    expect_true(absorbed[["null_diagnostics"]][["fit"]][["available"]])
    expect_equal(absorbed[["null_diagnostics"]][["weighted"]][["logLik"]],
      as.numeric(logLik(fit, type = "weighted")), tolerance = 1e-8)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})
