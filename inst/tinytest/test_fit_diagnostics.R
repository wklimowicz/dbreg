set.seed(830)
dat = data.frame(x = rnorm(300), z = rnorm(300), f = sample(1:9, 300, TRUE),
  g = sample(1:6, 300, TRUE), cl = sample(1:15, 300, TRUE), w = runif(300, 0.2, 2))
dat[["d"]] = dat[["z"]] + rnorm(300)
dat[["y"]] = dat[["x"]] + dat[["d"]] + dat[["f"]] + rnorm(300)
dat[["x"]][1:3] = NA_real_
dat[["cl"]][4:5] = NA_integer_
dat[["w"]][6] = 0
before = options()[c("dbreg.ap_tol", "dbreg.ap_max_iter")]
for (iv in c(FALSE, TRUE)) {
  fit = if (iv) {
    dbivreg(y ~ x | f + g | d ~ z, data = dat, vcov = ~cl, weights = "w",
      diagnostics = TRUE, fe_tol = 1e-7, fe_max_iter = 200)
  } else {
    dbreg(y ~ x + d | f + g, data = dat, vcov = ~cl, weights = "w",
      strategy = "demean", diagnostics = TRUE, fe_tol = 1e-7, fe_max_iter = 200)
  }
  expect_equal(fit[["nobs_input"]], 300)
  expect_equal(fit[["nobs_dropped"]], 6)
  expect_equal(fit[["sample_diagnostics"]][["n_nonpositive_weight"]], 1)
  detail = fit[["sample_diagnostics"]][["detail"]]
  expect_equal(detail[["n_missing"]][detail[["variable"]] == "x"], 3)
  expect_equal(generics::glance(fit)[["vcov.type"]], "cluster")
}
fit = dbreg(y ~ x | f + g, data = dat, strategy = "demean", fe_tol = 1e-7, fe_max_iter = 200)
expect_equal(fit[["ap_tolerance"]], 1e-7)
expect_equal(fit[["ap_max_iter"]], 200L)
expect_true(is.na(fit[["nobs_input"]]))
expect_equal(options()[names(before)], before)
for (value in list(0, -1, NA_real_, Inf, c(1e-7, 1e-8))) {
  expect_error(dbreg(y ~ x, data = dat, fe_tol = value), "fe_tol")
}
for (value in list(0, 1.5, 2^32, NA_real_)) {
  expect_error(dbreg(y ~ x, data = dat, fe_max_iter = value), "fe_max_iter")
}
expect_error(dbreg(y ~ x, data = dat, diagnostics = NA), "diagnostics")
dat[["x_copy"]] = dat[["x"]]
fit = dbreg(y ~ x + x_copy | f, data = dat, strategy = "demean")
expect_equal(sum(is.na(coef(fit, collin = TRUE))), 1L)
tidied = generics::tidy(fit, collin = TRUE, conf.int = TRUE)
expect_equal(sum(is.na(tidied[["estimate"]])), 1L)
expect_equal(sum(is.na(tidied[["conf.low"]])), 1L)
