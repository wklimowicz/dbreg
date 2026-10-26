local({
  set.seed(1032)
  n = 300L
  d = data.frame(x = sample(0:3, n, TRUE), g = sample(1:5, n, TRUE),
    h = sample(1:3, n, TRUE), cl = sample(1:9, n, TRUE),
    w = stats::runif(n, 0.3, 1.2))
  d[["o"]] = 0.1 * (d[["x"]] %% 2L)
  for (family in list(stats::gaussian(), stats::poisson(), stats::binomial(),
    stats::binomial("probit"), stats::binomial("cloglog"))) {
    eta = 0.2 * d[["x"]] + d[["o"]]
    d[["y"]] = switch(family[["family"]],
      gaussian = eta + stats::rnorm(n),
      poisson = stats::rpois(n, family[["linkinv"]](eta)),
      binomial = stats::rbinom(n, 1L, family[["linkinv"]](eta)))
    for (kind in list("iid", "hc1", ~cl)) {
      fit = dbglm(y ~ x | g + h, data = d, family = family,
        weights = ~w, offset = ~o, vcov = kind, strategy = "compress",
        control = dbglm_control(projection_method = "sql"), fixef_rm = "none")
      reference = fixest::feglm(y ~ x | g + h, data = d, family = family,
        weights = ~w, offset = ~o, vcov = kind, fixef.rm = "none",
        fixef.tol = 1e-10, glm.tol = 1e-12, glm.iter = 100L)
      expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)), tolerance = 1e-6)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-5)
      expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(reference)), tolerance = 1e-6)
      expect_true(fit[["n_cells"]] < n)
    }
  }
})
