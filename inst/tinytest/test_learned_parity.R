# Learned fixed effects: new-data predictions and fixef() match fixest.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(20261117)
  n = 800
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = sample(30, n, TRUE),
    h = sample(10, n, TRUE), w = runif(n, .5, 2))
  eta = .3*d$x - .2*d$z + rnorm(30, 0, .3)[d$g] + rnorm(10, 0, .3)[d$h]
  d$y = eta + rnorm(n); d$yp = rpois(n, exp(eta))
  nd = d[sample(n, 200), ]
  for (rhs in c("x + z | g", "x + z | g + h", "x | g[z]")) for (s in c("auto", "compress")) {
    f = as.formula(paste("y ~", rhs)); info = paste(rhs, s)
    ref = fixest::feols(f, d, weights = ~w, fixef.tol = 1e-11, fixef.rm = "none")
    fit = dbreg(f, data = d, weights = "w", strategy = s, retain = "prediction")
    expect_equal(unname(predict(fit, nd)), unname(predict(ref, nd)), tolerance = 1e-8, info = info)
    if (rhs == "x + z | g") expect_equal(fixest::fixef(fit)[["g"]][names(fixest::fixef(ref)$g)],
      fixest::fixef(ref)$g, tolerance = 1e-8, info = info)
    close(fit)
  }
})
