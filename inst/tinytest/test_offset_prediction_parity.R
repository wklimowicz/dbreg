# Offsets and demeaned-FE predictions match fixest / lm.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(20261117)
  n = 800
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(150, n, TRUE), h = sample(10, n, TRUE),
    k = sample(5, n, TRUE), w = runif(n, .5, 2), o = runif(n, -.5, .5))
  d$y = .5*d$x - .2*d$z + rnorm(150)[d$g] + rnorm(10)[d$h] + rnorm(5)[d$k] + d$o + rnorm(n)
  for (f in list(y ~ x + z, y ~ x + z | g, y ~ x + z | g + h)) for (s in c("auto", "compress")) {
    ref = fixest::feols(f, d, offset = ~o, weights = ~w, fixef.rm = "none")
    fit = dbreg(f, data = d, offset = "o", weights = "w", strategy = s)
    expect_equal(coef(fit)[c("x", "z")], coef(ref)[c("x", "z")], tolerance = 1e-9, info = paste(deparse(f), s))
    expect_equal(sqrt(diag(as.matrix(vcov(fit))))[c("x", "z")], fixest::se(ref)[c("x", "z")],
      tolerance = 1e-8, info = paste(deparse(f), s))
  }
  # Demeaned predictions with two or three weighted families recover the full
  # least-squares fit (many small groups make projection converge slowly).
  for (rhs in c("g + h", "g + h + k")) {
    fit = dbreg(as.formula(paste("y ~ x + z |", rhs)), data = d, weights = "w", strategy = "demean")
    ref = lm(as.formula(paste("y ~ x + z +", gsub("(\\w)", "factor(\\1)", rhs))), d, weights = w)
    expect_equal(predict(fit, d), unname(fitted(ref)), tolerance = 1e-9, info = rhs)
  }
})
