if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(114)
  d = expand.grid(x = -1:1, g = factor(1:3), h = factor(1:4), rep = 1:4)
  d$w = runif(nrow(d), 0.1, 8)
  d$o = d$rep / 5
  d$z = factor(d$rep %% 3)
  d$y = d$x + d$o + as.integer(d$g) / 2 + rnorm(nrow(d))
  d$x[5L] = NA_real_
  ref = lm(y ~ x * z + g + h, data = d, weights = w, offset = o)
  pred = d[c(9, 17, 5, 30, 100), c("x", "z", "g", "h", "o")]
  for (strategy in c("demean", "compress", "mundlak")) {
    fit = dbreg(y ~ x * z | g + h, data = d, weights = ~w, offset = ~o,
      retain = "all", strategy = strategy)
    expected = fitted(fit)
    expect_equal(predict(fit, d)[!is.na(d$x)], expected, tolerance = 1e-10)
    if (strategy != "mundlak") {
      expect_equal(predict(fit, pred), unname(predict(ref, pred)), tolerance = 1e-10)
    }
    before = predict(fit, pred)
    close(fit)
    expect_equal(predict(unserialize(serialize(fit, NULL)), pred), before)
    expect_equal(predict(fit, pred[1, , drop = FALSE]), before[1])
    unseen = pred
    unseen$g = as.character(unseen$g)
    unseen$g[1L] = "unknown"
    expect_true(is.na(predict(fit, unseen)[1L]))
    expect_error(predict(fit, pred, interval = "confidence"), "interval")
  }

  d = expand.grid(g = factor(1:3), h = factor(1:2), trend = -1:2, x = -1:1, rep = 1:3)
  d$w = runif(nrow(d), 1, 3)
  d$y = 3 * d$x + 0.2 * as.integer(d$g) * d$trend + rnorm(nrow(d))
  for (strategy in c("demean", "compress")) {
    fit = dbreg(y ~ x | (g^h)[trend], data = d, weights = ~w, strategy = strategy, retain = "all")
    ref = lm(y ~ x + interaction(g, h) * trend, data = d, weights = w)
    expected = unname(predict(ref, d))
    close(fit)
    expect_equal(predict(fit, d[, c("g", "h", "x", "trend")]), expected, tolerance = 1e-9)
  }

  d = data.frame(g = rep(c("a", "b"), each = 6), h = rep(c("c", "d"), each = 6),
    x = rep(-1:1, 4), y = rnorm(12))
  d$copy = d$x
  fit = dbreg(y ~ x + copy | g + h, data = d, strategy = "demean", retain = "prediction")
  close(fit)
  pred = d[1:3, ]
  pred$h[2] = "d"
  pred$copy[3] = 4
  expect_equal(is.na(predict(fit, pred)), c(FALSE, TRUE, TRUE))

  d$flag = rep(c(FALSE, TRUE), 6)
  flag = dbreg(y ~ x + flag | g, data = d, strategy = "compress", retain = "prediction")
  close(flag)
  expect_equal(predict(flag, d), unname(fitted(lm(y ~ x + flag + g, data = d))), tolerance = 1e-10)

  previous = options(dbreg.prediction_max_entries = 1)
  limited = dbreg(y ~ x | g, data = d, strategy = "demean", retain = "prediction")
  options(previous)
  expect_identical(limited[["prediction_diagnostics"]][["reason"]], "system_limit")
  expect_error(predict(limited, d), "no portable learned-effect state")
  close(limited)
}
