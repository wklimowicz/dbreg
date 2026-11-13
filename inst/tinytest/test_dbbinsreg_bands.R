# Independent bin-wise Gaussian regressions validate the plotted covariance.
local({
  set.seed(731)
  n = 500L
  B = 5L
  nsims = 97L
  d = data.frame(x = (seq_len(n)-.5)/n, g = rep(1:10, length.out = n))
  d[["y"]] = 2+.4*d[["x"]]+rnorm(n, sd = .2+d[["x"]])
  breaks = seq(0, 1, length.out = B+1L)
  bin = cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE)
  indicators = stats::model.matrix(~bin-1)
  left = as.numeric(tapply(d[["x"]], bin, min))
  right = as.numeric(tapply(d[["x"]], bin, max))
  u = (d[["x"]]-left[as.integer(bin)])/(right-left)[as.integer(bin)]
  u_mean = as.numeric(tapply(u, bin, mean))
  for (degree in 0:1) {
    X = if (degree == 0L) indicators else cbind(indicators, indicators*u)
    Z = if (degree == 0L) diag(B) else cbind(diag(B), diag(u_mean))
    reference = stats::lm.fit(X, d[["y"]])
    bread = solve(crossprod(X))
    scores = X*reference[["residuals"]]
    V = n/(n-ncol(X))*(bread %*% crossprod(scores) %*% bread)
    covariance = Z %*% V %*% t(Z)
    standard_errors = sqrt(diag(covariance))
    fitted = as.numeric(Z %*% reference[["coefficients"]])
    set.seed(719)
    actual = dbbinsreg(y ~ x, data = d, points = c(degree, 0), binspos = breaks,
      vcov = "hc1", cb = TRUE, nsims = nsims, plot = FALSE)
    set.seed(719)
    draws = MASS::mvrnorm(nsims, mu = rep(0, B), Sigma = covariance)
    critical = as.numeric(stats::quantile(apply(abs(sweep(draws, 2L, standard_errors, "/")), 1L, max), .95))
    expect_equal(actual[["points"]][["fit"]], fitted, tolerance = 1e-8)
    expect_equal(actual[["points"]][["se"]], standard_errors, tolerance = 1e-8)
    expect_equal(actual[["points"]][["cb_lwr"]], fitted-critical*standard_errors, tolerance = 1e-7)
    expect_equal(actual[["points"]][["cb_upr"]], fitted+critical*standard_errors, tolerance = 1e-7)
    set.seed(719)
    bands = dbbinsreg(y ~ x, data = d, points = c(degree, 0), binspos = breaks,
      vcov = "hc1", ci = FALSE, cb = TRUE, nsims = nsims, plot = FALSE)
    expect_false("se" %in% names(bands[["points"]]))
    expect_equal(bands[["points"]][["cb_lwr"]], actual[["points"]][["cb_lwr"]], tolerance = 1e-8)
  }
  clustered = dbbinsreg(y ~ x, data = d, binspos = breaks, vcov = ~g, plot = FALSE)
  expect_equal(clustered[["points"]][["upr"]]-clustered[["points"]][["fit"]],
    stats::qt(.975, 9)*clustered[["points"]][["se"]], tolerance = 1e-8)
  for (count in c(0, 1, NA_real_)) expect_error(dbbinsreg(y ~ x, data = d, cb = TRUE, nsims = count, plot = FALSE), pattern = "nsims")
})
