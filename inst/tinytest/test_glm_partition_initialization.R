# Equivalent constant-spanning encodings need the same offset-aware start.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(y = rep(c(1, 2, 3, 2, 3, 4), 3),
    g = rep(c("a", "a", "a", "b", "b", "b"), 3))
  d[["z"]] = rep(c(-1, 0, 1), 6)
  for (location in c(-30, 30)) {
    d[["o"]] = location + .1 * d[["z"]]
    for (fml in list(y ~ 0 + g, y ~ 0 + g + z)) {
      reference = stats::glm(fml, data = d, offset = o, family = poisson(),
        control = stats::glm.control(epsilon = 1e-13))
      for (strategy in c("moments", "compress")) {
        fit = dbglm(fml, data = d, conn = conn, family = poisson(), offset = "o",
          strategy = strategy, retain = "inference")
        expect_true(fit[["converged"]] && fit[["finite_mle"]])
        expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
        expect_equal(predict(fit, d), unname(predict(reference)), tolerance = 1e-8)
        expect_equal(deviance(fit), deviance(reference), tolerance = 1e-9)
        X = stats::model.matrix(reference)
        expected = solve(crossprod(X, X * stats::fitted(reference))) *
          (nrow(d) - 1) / (nrow(d) - ncol(X))
        expect_equal(as.numeric(vcov(fit)), as.numeric(expected), tolerance = 1e-8)
        close(fit)
      }
    }
  }
  # Recovery must target the same constant direction as default initialization.
  d[["o"]] = 200
  d[["y"]] = rep(c(0, 1, 1, 1, 0, 0), 3)
  d[["y"]][7:12] = 1 - d[["y"]][7:12]
  reference = stats::glm(y ~ 0 + g + z, data = d, offset = o, family = binomial(),
    control = stats::glm.control(epsilon = 1e-13))
  for (strategy in c("moments", "compress")) {
    fit = dbglm(y ~ 0 + g + z, data = d, conn = conn, offset = "o", family = binomial(),
      strategy = strategy, start = c(ga = 1000, gb = 1000, z = 1000))
    expect_true(fit[["converged"]] && fit[["step_halvings"]] > 0)
    expect_equal(coef(fit), coef(reference), tolerance = 1e-7)
    expect_equal(predict(fit, d), unname(predict(reference)), tolerance = 1e-7)
    close(fit)
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
