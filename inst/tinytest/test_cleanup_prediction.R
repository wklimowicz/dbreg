local({
  set.seed(917)
  d = data.frame(period = rep(1:4, 50), treat = rnorm(200), fe1 = rep(1:10, each = 20))
  d[["y"]] = rnorm(200) + d[["treat"]] + d[["fe1"]] / 10
  fml = y ~ i(period, treat, ref = 1) | fe1
  fit = dbreg::dbreg(fml, data = d, retain = "prediction")
  reference = fixest::feols(fml, data = d, fixef.rm = "none", notes = FALSE)
  expect_equal(stats::predict(fit, newdata = d),
    stats::predict(reference, newdata = d), tolerance = 1e-12)
  expect_equal(stats::predict(fit, newdata = d[200:1, ]),
    stats::predict(reference, newdata = d[200:1, ]), tolerance = 1e-12)
})
