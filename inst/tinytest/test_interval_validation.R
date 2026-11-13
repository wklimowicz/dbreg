# Invalid confidence levels fail early before reading prediction inputs.
local({
  set.seed(531)
  d = data.frame(x = rnorm(120), g = factor(rep(1:6, each = 20)))
  d[["y"]] = 0.5 * d[["x"]] + as.numeric(d[["g"]]) + rnorm(120)

  plain = dbreg(y ~ x, data = d, strategy = "moments")
  bad_levels = list(0, 1, 2, -0.1, NA_real_, Inf, -Inf, c(0.9, 0.95), "95%", numeric(0))
  for (value in bad_levels) {
    expect_error(confint(plain, level = value), "level")
    for (interval in c("confidence", "prediction")) {
      expect_error(predict(plain, newdata = NULL, interval = interval, level = value), "level")
      expect_error(predict(plain, data.frame(), interval = interval, level = value), "level")
    }
  }

  # Valid non-default intervals match lm() for plain and explicit-FE models.
  dense_plain = lm(y ~ x, data = d)
  expect_equal(unname(confint(plain, level = 0.9)), unname(confint(dense_plain, level = 0.9)), tolerance = 1e-7)
  for (interval in c("confidence", "prediction")) {
    expected = predict(dense_plain, d, interval = interval, level = 0.9)
    actual = predict(plain, d, interval = interval, level = 0.9)
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = 1e-7)
  }

  fe_mod = dbreg(y ~ x | g, data = d, strategy = "compress")
  dense_fe = lm(y ~ x + g, data = d)
  expect_equal(
    unname(confint(fe_mod, level = 0.9)),
    unname(confint(dense_fe, parm = "x", level = 0.9)),
    tolerance = 1e-7
  )
  for (interval in c("confidence", "prediction")) {
    expected = predict(dense_fe, d, interval = interval, level = 0.9)
    actual = predict(fe_mod, d, interval = interval, level = 0.9)
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = 1e-7)
  }

  # Legacy `fes` alias continues to work with a deprecation warning.
  expect_warning({ with_fe = coef(fe_mod, fes = TRUE) }, "superseded")
  expect_equal(with_fe, coef(fe_mod, fe = TRUE))
  expect_warning({ ci = confint(fe_mod, fes = TRUE) }, "superseded")
  expect_equal(ci, confint(fe_mod, fe = TRUE))
  expect_warning({ full = vcov(fe_mod, fes = TRUE) }, "superseded")
  expect_equal(full, vcov(fe_mod, fe = TRUE))
})
