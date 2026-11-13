local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11164)
  d = expand.grid(g = 1:5, x = -1:1, rep = 1:24)
  d$w = runif(nrow(d), .4, 2)
  d$o = sin(d$rep) / 10
  d$cluster = d$rep
  cases = list(list(gaussian(), "fisher"), list(poisson(), "fisher"),
    list(binomial(), "fisher"), list(binomial("probit"), "fisher"),
    list(binomial("cloglog"), "fisher"), list(binomial("cloglog"), "newton"))
  for (case in cases) {
    family = case[[1L]]
    eta = .1 + .35 * d$x + cos(d$g) / 3 + d$o
    mu = family$linkinv(eta)
    d$y = switch(family$family, gaussian = mu + rnorm(nrow(d)),
      poisson = rpois(nrow(d), mu), binomial = rbinom(nrow(d), 1, mu))
    dense = suppressWarnings(glm(y ~ x + factor(g), data = d, family = family,
      weights = w, offset = o, control = glm.control(epsilon = 1e-15, maxit = 1000L)))
    for (strategy in c("moments", "compress")) {
      fit = dbglm(y ~ x | g, data = d, conn = conn, family = family,
        weights = ~w, offset = ~o, strategy = strategy, retain = "all",
        control = dbglm_control(projection_method = "sql", irls = case[[2L]], score_tol = 1e-9))
      expect_identical(fit$projection_method, "sql_window")
      expect_true(fit$converged && max(fit$score_max, fit$fe_score_max) <= 1e-9)
      expect_equal(unname(coef(fit)), unname(coef(dense)["x"]), tolerance = 1e-8)
      expect_equal(unname(predict(fit, d, type = "link")), unname(predict(dense, type = "link")), tolerance = 1e-8)
      # Base glm supplies the unadjusted expected information (and Gaussian
      # dispersion); the package's separate requested SSC is applied explicitly.
      K = qr(model.matrix(dense))$rank
      correction = (nrow(d) - 1) / (nrow(d) - K)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(dense)["x", "x"]) * correction,
        tolerance = 1e-8)
      X = model.matrix(dense)
      eta_dense = predict(dense, type = "link")
      W = d$w * family$mu.eta(eta_dense)^2 / family$variance(fitted(dense))
      B = solve(crossprod(X, X * W))
      scores = X * (d$w * (d$y - fitted(dense)) * family$mu.eta(eta_dense) / family$variance(fitted(dense)))
      clustered = rowsum(scores, d$cluster)
      expected = B %*% crossprod(clustered) %*% B
      observed = vcov(fit, cluster = ~cluster, ssc = list(K.adj = FALSE, G.adj = FALSE))
      expect_equal(as.numeric(observed), as.numeric(expected["x", "x"]), tolerance = 1e-8)
      close(fit)
      expect_identical(DBI::dbListTables(conn), character())
    }
  }
  d$warm = unname(predict(dense, type = "link"))
  fit = dbglm(y ~ x | g, data = d, conn = conn, family = family,
    weights = ~w, offset = ~o, etastart = ~warm,
    control = dbglm_control(projection_method = "sql", compute_inference = FALSE))
  expect_true(fit$converged)
  expect_equal(fit$iterations, 1L)
  expect_true(all(is.na(vcov(fit))))
  close(fit)
  expect_identical(DBI::dbListTables(conn), character())
})
