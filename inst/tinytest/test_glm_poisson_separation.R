# A disappearing Poisson score is not a finite solution when a coefficient
# can decrease zero-count means indefinitely without changing positive rows.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(x = rep(0:1, each = 20),
    y = c(rep(c(1, 2, 3, 2), 5), rep(0, 20)))
  d[["u"]] = rep(c(-1, 1), 20)
  d[["v"]] = d[["x"]] - d[["u"]]
  initial = sort(DBI::dbListTables(conn))
  for (strategy in c("moments", "compress")) {
    for (fml in list(y ~ x, y ~ u + v)) {
      expect_error(dbglm(fml, data = d, conn = conn, family = poisson(), strategy = strategy), "finite|boundary|separation")
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
    expect_error(dbglm(y ~ x, data = d, conn = conn, family = poisson(),
      strategy = strategy, start = c("(Intercept)" = log(2), x = -30)), "finite|boundary|separation")
    expect_error(dbglm(y ~ x, data = d, conn = conn, family = poisson(),
      strategy = strategy, separation = "prediction"), "binomial")
    zero = data.frame(y = rep(0, 40), x = rep(1:2, each = 20))
    expect_error(dbglm(y ~ 0 + x, data = zero, conn = conn, family = poisson(), strategy = strategy), "finite|boundary|separation")
    # A small but positive weighted group count still has a finite MLE, even
    # when its information decreases sharply from the initial design.
    finite = d
    finite[["y"]][21L] = 1
    finite[["w"]] = 1
    finite[["w"]][21L] = 1e-12
    mean_zero = finite[["w"]][21L] / sum(finite[["w"]][finite[["x"]] == 1])
    expected = c(log(2), log(mean_zero) - log(2))
    fit = dbglm(y ~ x, data = finite, conn = conn, weights = "w", family = poisson(), strategy = strategy)
    expect_true(fit[["finite_mle"]] && fit[["converged"]])
    # The information condition number is about 1e14 in this finite control.
    expect_equal(unname(coef(fit)), expected, tolerance = 1e-4)
    expect_equal(predict(fit, finite, type = "response"), ifelse(finite[["x"]] == 0, 2, mean_zero), tolerance = 1e-8)
    close(fit)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})
