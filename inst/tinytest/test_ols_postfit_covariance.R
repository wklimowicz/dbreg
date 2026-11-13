local({
  set.seed(11156)
  n = 300L
  d = data.frame(g = rep(1:6, each = 50), h = rep(1:10, length.out = n),
    cl = sample(1:15, n, TRUE), x = sample(-2:2, n, TRUE), w = runif(n, .1, 3))
  d[["y"]] = .4 * d[["x"]] + sin(d[["g"]]) + rnorm(n)
  d[["bad"]] = d[["cl"]]
  d[["bad"]][1L] = NA
  ref = lm(y ~ x + factor(g), data = d, weights = w)
  X = model.matrix(ref)
  bread = solve(crossprod(X, X * d[["w"]]))
  score = X * (d[["w"]] * residuals(ref))
  p = ncol(X)
  cross = function(g) crossprod(rowsum(score, g))
  G = length(unique(d[["cl"]]))
  H = length(unique(d[["h"]]))
  adjust = min(G, H) / (min(G, H) - 1) * (n - 1) / (n - p)
  two = bread %*% (cross(d[["cl"]]) + cross(d[["h"]]) -
    cross(interaction(d[["cl"]], d[["h"]]))) %*% bread * adjust
  single = bread %*% cross(d[["cl"]]) %*% bread * G / (G - 1) * (n - 1) / (n - p)
  hc1 = bread %*% crossprod(score) %*% bread * n / (n - p)
  for (strategy in c("demean", "compress")) {
    fit = dbreg(y ~ x | g, data = d, weights = ~w, retain = "inference",
      strategy = strategy, vcov = "iid")
    expect_equal(vcov(fit, vcov = "hc1")["x", "x"], hc1["x", "x"], tolerance = 1e-9)
    expect_equal(vcov(fit, vcov = ~cl)["x", "x"], single["x", "x"], tolerance = 1e-9)
    expect_equal(suppressWarnings(vcov(fit, vcov = ~cl + h, vcov_fix = FALSE))["x", "x"],
      two["x", "x"], tolerance = 1e-9)
    expect_equal(vcov(fit, vcov = "iid")["x", "x"], vcov(ref)["x", "x"], tolerance = 1e-9)
    expect_error(vcov(fit, vcov = ~bad), "missing or nonfinite")
    expect_equal(attr(vcov(fit, vcov = ~cl), "df.t"), G - 1L)
    close(fit)
    expect_error(vcov(fit, vcov = "hc1"), "no live database session")
    expect_true(is.finite(vcov(fit)["x", "x"]))
  }
})
