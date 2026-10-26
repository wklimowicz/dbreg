local({
  set.seed(1026)
  n = 180L
  d = data.frame(y = stats::rnorm(n), g = sample(1:6, n, TRUE),
    z = sample(0:2, n, TRUE), h = sample(1:9, n, TRUE),
    w = stats::runif(n, 0.3, 1.5))
  names = c("n", "sum_w", "sum_wy", "sum_wy_sq", "sum_w2", "sum_w2y",
    "sum_w2y_sq", "sum_w_gc", "sum_wy_gc", "mean_Y", "wts",
    "dbreg_cell_center__", "dbreg_cell_scale__")
  for (weighted in c(FALSE, TRUE)) {
    weights = if (weighted) ~w else NULL
    reference = fixest::feols(y ~ z | g, data = d, weights = weights,
      vcov = "hc1", fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    clustered = fixest::feols(y ~ z | g, data = d, weights = weights,
      vcov = ~h, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    for (name in names) {
      data = d
      data[[name]] = data[["z"]]
      formula = stats::as.formula(paste("y ~", name, "| g"))
      fit = dbreg(formula, data = data, weights = weights,
        strategy = "compress", vcov = "hc1", retain = "all")
      expect_identical(names(stats::coef(fit)), name)
      expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)), tolerance = 1e-7)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
      expect_equal(stats::predict(fit, data), as.numeric(stats::predict(reference)), tolerance = 1e-7)
      expect_equal(as.numeric(stats::vcov(fit, vcov = ~h)), as.numeric(stats::vcov(clustered)), tolerance = 1e-7)
      close(fit)
    }
  }
})
