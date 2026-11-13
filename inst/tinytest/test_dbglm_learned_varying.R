local({
  set.seed(11153)
  d = expand.grid(g = 1:3, h = 1:2, x = c(-1, 1), z = c(-1, 1), rep = 1:7)
  d[["trend"]] = 1e12 + d[["z"]]
  d[["constant"]] = 5
  d[["w"]] = runif(nrow(d), .3, 3)
  d[["y"]] = .3 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) * d[["z"]] + rnorm(nrow(d))
  fit = dbglm(y ~ x | (g^h)[trend, constant], data = d, weights = ~w,
    retain = "prediction", strategy = "compress")
  expect_identical(fit[["prediction_diagnostics"]][["status"]], "completed")
  expect_error(fitted(fit), "inference")
  ref = lm(y ~ x + interaction(g, h)*z, data = d, weights = w)
  expect_equal(predict(fit, d), unname(fitted(ref)), tolerance = 1e-10)
  new = d[rep(1L, 4L), ]
  new[["z"]][1L] = .5
  new[["trend"]][1L] = 1e12 + .5
  new[["constant"]][2L] = 6
  new[["g"]][3L] = 4
  new[["trend"]][4L] = NA_real_
  close(fit)
  expect_identical(is.na(predict(fit, new)), c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(unname(predict(fit, new[1L, ])), unname(predict(ref, new[1L, ])), tolerance = 1e-10)
  expect_equal(predict(unserialize(serialize(fit, NULL)), d),
    unname(fitted(ref)), tolerance = 1e-10)
})
