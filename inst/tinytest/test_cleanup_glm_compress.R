local({
  set.seed(1012)
  n = 2400L
  d = data.frame(x = sample(0:3, n, TRUE), g = sample(1:8, n, TRUE),
    h = sample(1:4, n, TRUE), w = sample(c(0.2, 0.5, 1, 1.5, 2.5), n, TRUE))
  d[["o"]] = 0.1 * (d[["x"]] %% 2)
  eta = 0.1 + 0.2 * d[["x"]] + sin(d[["g"]]) * 0.15 + cos(d[["h"]]) * 0.1 + d[["o"]]
  for (family in list(stats::gaussian(), stats::poisson(), stats::binomial("logit"),
    stats::binomial("probit"), stats::binomial("cloglog"))) {
    mu = family[["linkinv"]](eta)
    d[["y"]] = switch(family[["family"]], gaussian = mu + stats::rnorm(n),
      poisson = stats::rpois(n, mu), binomial = stats::rbinom(n, 1L, mu))
    for (fml in list(y ~ x, y ~ x | g + h)) for (kind in list("iid", "hc1", ~g)) {
      ordinary = dbglm(fml, data = d, weights = ~w, offset = ~o, family = family,
        vcov = kind, strategy = "moments", fixef_rm = "none")
      compressed = dbglm(fml, data = d, weights = ~w, offset = ~o, family = family,
        vcov = kind, strategy = "compress", fixef_rm = "none")
      expect_equal(stats::coef(compressed), stats::coef(ordinary), tolerance = 1e-7)
      expect_equal(as.numeric(stats::vcov(compressed)), as.numeric(stats::vcov(ordinary)), tolerance = 1e-7)
      expect_equal(as.numeric(stats::logLik(compressed)), as.numeric(stats::logLik(ordinary)), tolerance = 1e-8)
      expect_equal(as.numeric(stats::logLik(compressed, type = "weighted")),
        as.numeric(stats::logLik(ordinary, type = "weighted")), tolerance = 1e-8)
      expect_equal(compressed[["deviance"]], ordinary[["deviance"]], tolerance = 1e-8)
      expect_identical(stats::nobs(compressed), stats::nobs(ordinary))
      expect_identical(compressed[["n_params_model"]], ordinary[["n_params_model"]])
      expect_identical(compressed[["strategy"]], "compress")
      expect_true(compressed[["n_cells"]] < n / 4)
    }
    ordinary = dbglm(y ~ x | g + h, data = d, weights = ~w, offset = ~o, family = family,
      strategy = "moments", retain = "prediction", fixef_rm = "none")
    compressed = dbglm(y ~ x | g + h, data = d, weights = ~w, offset = ~o, family = family,
      strategy = "compress", retain = "prediction", fixef_rm = "none")
    expect_equal(stats::predict(compressed, d), stats::predict(ordinary, d), tolerance = 1e-7)
    expect_equal(stats::predict(compressed, d, type = "response"),
      stats::predict(ordinary, d, type = "response"), tolerance = 1e-7)
    retained_data = d
    retained_data[["cl"]] = sample(1:9, n, TRUE)
    retained_data[["y"]][1L] = NA_real_
    retained_data[["w"]][2L] = 0
    for (fml in list(y ~ x, y ~ x | g + h)) for (retain in c("inference", "all")) {
      ordinary = dbglm(fml, data = retained_data, weights = ~w, offset = ~o, family = family,
        strategy = "moments", retain = retain, fixef_rm = "none")
      compressed = dbglm(fml, data = retained_data, weights = ~w, offset = ~o, family = family,
        strategy = "compress", retain = retain, fixef_rm = "none")
      for (kind in list("iid", "hc1", ~cl)) {
        expect_equal(as.numeric(stats::vcov(compressed, vcov = kind)),
          as.numeric(stats::vcov(ordinary, vcov = kind)), tolerance = 1e-7)
      }
      expect_equal(stats::predict(compressed), stats::predict(ordinary), tolerance = 1e-7)
      expect_equal(stats::predict(compressed, sample = "original"),
        stats::predict(ordinary, sample = "original"), tolerance = 1e-7)
      expect_equal(stats::residuals(compressed, na.rm = FALSE),
        stats::residuals(ordinary, na.rm = FALSE), tolerance = 1e-7)
      expect_identical(stats::case.names(compressed), stats::case.names(ordinary))
      session = dbreg:::fitted_resource_connection(compressed[["resource_state"]])
      expect_equal(length(DBI::dbListTables(session)), 1L)
      close(ordinary)
      close(compressed)
    }
  }
})
