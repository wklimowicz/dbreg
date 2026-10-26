local({
  set.seed(1013)
  n = 480L
  d = data.frame(x = sample(0:3, n, TRUE), g = sample(1:6, n, TRUE), h = sample(1:4, n, TRUE))
  weights = sample(c(0.2, 0.5, 1, 1.5, 2.5), n, TRUE)
  d[["o"]] = 0.1 * (d[["x"]] %% 2)
  eta = 0.15 * d[["x"]] + 0.1 * sin(d[["g"]]) + 0.05 * cos(d[["h"]]) + d[["o"]]
  for (family in list(stats::gaussian(), stats::poisson(), stats::binomial("logit"),
    stats::binomial("probit"), stats::binomial("cloglog"))) {
    message("[GLM weight scale] ", family[["family"]], "/", family[["link"]])
    mu = family[["linkinv"]](eta)
    d[["y"]] = switch(family[["family"]], gaussian = mu + stats::rnorm(n),
      poisson = stats::rpois(n, mu), binomial = stats::rbinom(n, 1L, mu))
    for (fml in list(y ~ x, y ~ x | g + h)) for (kind in list("iid", "hc1", ~g)) {
      d[["w"]] = weights
      reference = dbglm(fml, data = d, weights = ~w, offset = ~o, family = family,
        vcov = kind, strategy = "moments", fixef_rm = "none")
      for (scale in c(1e-200, 1e200)) for (strategy in c("moments", "compress")) {
        d[["w"]] = weights * scale
        fit = dbglm(fml, data = d, weights = ~w, offset = ~o, family = family,
          vcov = kind, strategy = strategy, fixef_rm = "none")
        expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-7)
        V = stats::vcov(fit)
        if (family[["family"]] != "gaussian" && identical(kind, "iid")) V = V * scale
        expect_equal(as.numeric(V), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
        expect_equal(fit[["deviance"]] / scale, reference[["deviance"]], tolerance = 1e-8)
        likelihood = as.numeric(stats::logLik(fit, type = "weighted"))
        if (family[["family"]] != "gaussian") likelihood = likelihood / scale
        expect_equal(likelihood, as.numeric(stats::logLik(reference, type = "weighted")), tolerance = 1e-8)
        if (family[["family"]] == "binomial") {
          if (scale < 1) expect_equal(as.numeric(stats::logLik(fit)), 0) else
            expect_equal(as.numeric(stats::logLik(fit)) / scale,
              as.numeric(stats::logLik(reference, type = "weighted")), tolerance = 1e-8)
        }
      }
    }
  }
})
