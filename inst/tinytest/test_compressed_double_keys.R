local({
  set.seed(11137)
  n = 180L
  group = rep(1:3, length.out = n)
  keys = 1 + (0:2) * .Machine[["double.eps"]]
  d = data.frame(g = keys[group], x = sample(c(-1, 1), n, TRUE),
    w = runif(n, .2, 2), cl = sample(1:10, n, TRUE))
  d[["y"]] = .6 * d[["x"]] + group + rnorm(n)
  X = model.matrix(~ x + factor(group), d)
  ref = lm.wfit(X, d[["y"]], d[["w"]])
  scores = rowsum(X * (d[["w"]] * ref[["residuals"]]), d[["cl"]])
  bread = solve(crossprod(X, d[["w"]] * X))
  G = nrow(scores)
  V = bread %*% crossprod(scores) %*% bread *
    (G / (G - 1)) * ((n - 1) / (n - ncol(X)))

  fit = dbreg(y ~ x | g, data = d, weights = ~w, vcov = ~cl,
    strategy = "compress")
  expect_equal(nobs(fit), n)
  expect_equal(fit[["nobs"]], 6)
  expect_equal(coef(fit)[["x"]], ref[["coefficients"]][["x"]], tolerance = 1e-10)
  expect_equal(vcov(fit)["x", "x"], V["x", "x"], tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, d)), as.numeric(ref[["fitted.values"]]),
    tolerance = 1e-10)
  cells = fit[["data"]]
  expect_identical(sort(unique(as.numeric(as.character(cells[["g"]])))), keys)
  expect_equal(as.numeric(predict(fit)), as.numeric(predict(fit, cells)),
    tolerance = 1e-10)
  restored = unserialize(serialize(fit, NULL))
  expect_equal(as.numeric(predict(restored, d)), as.numeric(ref[["fitted.values"]]),
    tolerance = 1e-10)
  unseen = d[1L, ]
  unseen[["g"]] = 1 + 3 * .Machine[["double.eps"]]
  expect_true(is.na(predict(restored, unseen)))
})
