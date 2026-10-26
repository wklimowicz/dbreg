local({
  set.seed(1005)
  n = 10000L
  source = data.frame(x = 1e9 + rnorm(n), z = rnorm(n), g = rep(1:50, length.out = n))
  source[["y"]] = 2e9 + 0.4 * (source[["x"]] - 1e9) - 0.2 * source[["z"]] + rnorm(n)
  reference_data = source
  reference_data[["x"]] = source[["x"]] - 1e9
  reference_data[["y"]] = source[["y"]] - 2e9
  for (kind in list("iid", "hc1", ~g)) {
    reference = fixest::feols(y ~ x + z, data = reference_data,
      vcov = if (identical(kind, "hc1")) "hetero" else kind)
    fit = dbreg::dbreg(y ~ x + z, data = source, strategy = "moments", vcov = kind)
    expect_equal(unname(stats::coef(fit)[-1L]), unname(stats::coef(reference)[-1L]), tolerance = 1e-6)
    expect_true(max(abs(stats::predict(fit, source) - (stats::predict(reference, reference_data) + 2e9))) < 1e-6)
    expect_equal(as.numeric(stats::vcov(fit)[-1L, -1L]), as.numeric(stats::vcov(reference)[-1L, -1L]), tolerance = 1e-6)
    expect_equal(unname(dbreg::gof(fit)[["rmse"]]), sqrt(mean(stats::residuals(reference)^2)), tolerance = 1e-6)
  }
  # Metadata-free NaN filtering must preserve dates and literal text "NaN".
  d = data.frame(x = rnorm(240), z = rnorm(240),
    text = rep(c("NaN", "ok"), 120), g = as.Date("2026-01-01") + rep(1:8, each = 30))
  d[["y"]] = 0.4 * d[["x"]] + (d[["text"]] == "ok") + rnorm(240)
  d[["x"]][1L] = NaN
  d[["z"]][2L] = NA_real_
  for (strategy in c("compress", "demean")) {
    fit = dbreg::dbreg(y ~ x + z + text | g, data = d, strategy = strategy)
    reference = fixest::feols(y ~ x + z + text | g, data = d, fixef.rm = "none")
    expect_equal(stats::nobs(fit), stats::nobs(reference))
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-7)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
  }
})
