local({
  set.seed(11138)
  N = 150L
  d = data.frame(n = rep(1:3, length.out = N),
    sum_w = sample(c(-1, 1), N, TRUE), w = runif(N, .1, 3),
    cl = sample(1:10, N, TRUE))
  d[["y"]] = .7 * d[["sum_w"]] + d[["n"]] + rnorm(N)
  X = model.matrix(~ sum_w + factor(n), d)
  ref = lm.wfit(X, d[["y"]], d[["w"]])
  bread = solve(crossprod(X, d[["w"]] * X))
  scores = rowsum(X * (d[["w"]] * ref[["residuals"]]), d[["cl"]])
  G = nrow(scores)
  V = bread %*% crossprod(scores) %*% bread *
    (G / (G - 1)) * ((N - 1) / (N - ncol(X)))

  fit = dbreg(y ~ sum_w | n, data = d, weights = ~w, vcov = ~cl,
    strategy = "compress")
  expect_identical(names(coef(fit)), "sum_w")
  expect_equal(coef(fit)[["sum_w"]], ref[["coefficients"]][["sum_w"]], tolerance = 1e-10)
  expect_equal(vcov(fit)["sum_w", "sum_w"], V["sum_w", "sum_w"], tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, d)), as.numeric(ref[["fitted.values"]]),
    tolerance = 1e-10)
  cells = dbreg(y ~ sum_w | n, data = d, weights = ~w, strategy = "compress",
    data_only = TRUE)
  map = attr(cells, "dbreg_source_columns")
  expect_identical(unname(map[c("n", "sum_w")]), c("n_source", "sum_w_source"))
  expect_equal(sum(cells[["n"]]), N)
  expect_equal(sum(cells[["sum_w"]]), sum(d[["w"]]))
  expect_equal(sum(cells[["sum_wy"]]), sum(d[["w"]] * d[["y"]]))
  expect_equal(as.numeric(predict(fit, cells)), as.numeric(predict(fit)),
    tolerance = 1e-10)
})
