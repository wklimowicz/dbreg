local({
  set.seed(11136)
  n = 400L
  d = data.frame(
    x = rnorm(n),
    g = sample(1:3, n, TRUE),
    w = runif(n, .2, 2),
    cluster = sample(1:12, n, TRUE)
  )
  d[["y"]] = .8 * d[["x"]] + sin(d[["g"]]) + rnorm(n)
  X = model.matrix(~ x + factor(g), d)
  ref = lm.wfit(X, d[["y"]], d[["w"]])
  bread = solve(crossprod(X, d[["w"]] * X))
  scores = rowsum(X * (d[["w"]] * ref[["residuals"]]), d[["cluster"]])
  G = nrow(scores)
  V = bread %*% crossprod(scores) %*% bread *
    (G / (G - 1)) * ((n - 1) / (n - ncol(X)))

  model = dbreg(y ~ x + g | g, data = d, weights = ~w, vcov = ~cluster,
    strategy = "compress")
  expect_identical(names(coef(model)), "x")
  expect_true("g" %in% model[["collin.var"]])
  expect_equal(coef(model)[["x"]], ref[["coefficients"]][["x"]], tolerance = 1e-10)
  expect_equal(vcov(model)["x", "x"], V["x", "x"], tolerance = 1e-10)
  expect_equal(as.numeric(predict(model, d)), as.numeric(ref[["fitted.values"]]),
    tolerance = 1e-10)
})
