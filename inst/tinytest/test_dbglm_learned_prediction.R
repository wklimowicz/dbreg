local({
  set.seed(11152)
  d = expand.grid(g = 1:4, h = 1:4, rep = 1:6)
  d = d[(d[["g"]] <= 2L) == (d[["h"]] <= 2L), ]
  d[["x"]] = rnorm(nrow(d))
  d[["z"]] = 2 * d[["x"]]
  d[["w"]] = runif(nrow(d), .2, 2)
  d[["o"]] = runif(nrow(d), -.2, .2)
  d[["y"]] = .35 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + d[["o"]] + rnorm(nrow(d))
  fit = dbglm(y ~ x + z | g + h, data = d, weights = ~w, offset = ~o,
    retain = "prediction", strategy = "moments")
  expect_identical(fit[["prediction_diagnostics"]][["status"]], "completed")
  expect_error(fitted(fit), "inference")
  ref = lm(y ~ x + factor(g) + factor(h), data = d, weights = w, offset = o)
  expect_equal(predict(fit, d), unname(fitted(ref)), tolerance = 1e-10)
  new = d[c(1L, 1L, 1L, 1L), ]
  new[["h"]][2L] = 3L
  new[["g"]][3L] = 20L
  new[["z"]][4L] = new[["z"]][4L] + 1
  expect_identical(is.na(predict(fit, new)), c(FALSE, TRUE, TRUE, TRUE))
  close(fit)
  expect_equal(predict(fit, d), unname(fitted(ref)), tolerance = 1e-10)
  restored = unserialize(serialize(fit, NULL))
  expect_equal(predict(restored, d), unname(fitted(ref)), tolerance = 1e-10)
  tuple = dbglm(y ~ x | g^h, data = d, weights = ~w, offset = ~o,
    retain = "all", control = dbglm_control(projection_method = "sql"))
  expect_equal(predict(tuple, d), fitted(tuple), tolerance = 1e-10)
  expect_true(is.na(predict(tuple, new[2L, ])))
  close(tuple)
  previous = options(dbreg.prediction_max_entries = 1)
  expect_warning(limited <- dbglm(y ~ x | g, data = d, retain = "all"), "could not be retained")
  options(previous)
  expect_identical(limited[["prediction_diagnostics"]][["reason"]], "system_limit")
  expect_equal(length(fitted(limited)), nrow(d))
  close(limited)
})
