# nobs() counts fitted source rows, not SQL result rows, cell counts or weights.
local({
  set.seed(901)
  d = data.frame(x = rep(-2:2, 24), g = rep(1:6, each = 20), w = rep(c(.5, 1, 2), 40))
  d[["z"]] = rnorm(120)
  d[["endog"]] = d[["z"]] + rnorm(120)
  d[["y"]] = .7 + .5*d[["x"]] + .2*d[["endog"]] + sin(d[["g"]]) + rnorm(120)
  d[["w"]][3] = 0; d[["x"]][4] = NA; d[["y"]][5] = NA; d[["w"]][6] = NA
  expected = sum(complete.cases(d) & d[["w"]] > 0, na.rm = TRUE)
  for (strategy in c("moments", "compress", "demean", "mundlak")) {
    fml = if (strategy %in% c("demean", "mundlak")) y ~ x | g else y ~ x
    fit = dbreg(fml, data = d, weights = "w", strategy = strategy)
    expect_identical(nobs(fit), as.numeric(expected))
    expect_equal(nobs(fit), fit[["nobs_orig"]])
    expect_true(nobs(fit) != sum(d[["w"]], na.rm = TRUE))
    if (strategy %in% c("moments", "compress")) expect_true(nobs(fit) > fit[["nobs"]])
  }
  for (strategy in c("moments", "compress", "demean")) {
    fml = if (strategy == "demean") y ~ x | g | endog ~ z else y ~ x | endog ~ z
    fit = dbivreg(fml, data = d, weights = "w", strategy = strategy)
    expect_identical(nobs(fit), as.numeric(expected))
    expect_equal(nobs(fit), fit[["nobs_orig"]])
  }
})
