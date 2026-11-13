# Retained residuals and post-fit covariance must consume the same stable SQL
# coordinates as the fit, even though reported coefficients use raw origins.
local({
  set.seed(867393)
  n = 160L
  d = data.frame(g = rep(1:8, each = 20), z = rnorm(n),
    x = 1e14 + runif(n, -1, 1), w = runif(n, .4, 2))
  d[["dx"]] = d[["x"]] - 1e14
  d[["y"]] = 2*d[["dx"]] + .4*d[["z"]] + sin(d[["g"]]) + rnorm(n, sd = .05)
  for (strategy in c("moments", "compress", "mundlak")) {
    reference = d
    if (strategy == "mundlak") {
      reference[["mx"]] = ave(seq_len(n), d[["g"]], FUN = function(i) weighted.mean(d[["dx"]][i], d[["w"]][i]))
      reference[["mz"]] = ave(seq_len(n), d[["g"]], FUN = function(i) weighted.mean(d[["z"]][i], d[["w"]][i]))
    }
    ref = lm(if (strategy == "moments") y ~ dx + z else if (strategy == "mundlak") y ~ dx + z + mx + mz else y ~ dx + z + factor(g),
      reference, weights = w)
    fit = dbreg(if (strategy == "moments") y ~ x + z else y ~ x + z | g,
      data = d, weights = "w", strategy = strategy, retain = "inference")
    expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(residuals(fit)), as.numeric(residuals(ref)), tolerance = 1e-8)
    expect_equal(predict(fit, d), as.numeric(fitted(ref)), tolerance = 1e-8)
    X = model.matrix(ref)
    bread = solve(crossprod(X, d[["w"]]*X))
    meat = crossprod(X*as.numeric(d[["w"]]*residuals(ref)))
    V = n/(n-ncol(X))*bread %*% meat %*% bread
    actual = vcov(fit, vcov = "hc1")
    positions = match(c("x", "z"), rownames(actual))
    expect_equal(as.numeric(actual[positions, positions]), as.numeric(V[c("dx", "z"), c("dx", "z")]), tolerance = 1e-8)
    close(fit)
  }
})
