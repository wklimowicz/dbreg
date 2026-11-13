# First-stage instrument strength is invariant to a common weight multiplier,
# including a restricted model with no explicit columns after FE absorption.
local({
  set.seed(998)
  n = 300L
  d = data.frame(z = rnorm(n), g = rep(1:15, each = 20))
  d[["d"]] = .4*d[["z"]] + rnorm(n)
  d[["y"]] = 2*d[["d"]] + rnorm(n)
  for (varying in c(FALSE, TRUE)) for (multiplier in c(.1, 1, 10, 1e6)) {
    d[["w"]] = multiplier * if (varying) seq(.5, 1.5, length.out = n) else 1
    for (absorbed in c(FALSE, TRUE)) {
      fit = if (absorbed) dbivreg(y ~ 1 | g | d ~ z, data = d, weights = "w") else
        dbivreg(y ~ 0 | d ~ z, data = d, weights = "w")
      # Direct weighted regressions independently supply both residual sums.
      full = if (absorbed) lm(d ~ z + factor(g), d, weights = w) else
        lm(d ~ 0 + z, d, weights = w)
      restricted = if (absorbed) sum(d[["w"]] * residuals(lm(d ~ factor(g), d, weights = w))^2) else
        sum(d[["w"]] * d[["d"]]^2)
      unrestricted = sum(d[["w"]] * residuals(full)^2)
      statistic = (restricted - unrestricted) / (unrestricted / df.residual(full))
      actual = fit[["diagnostics"]][["first_stage_f"]][["d"]]
      expect_equal(actual[["stat"]], statistic, tolerance = 1e-9)
      expect_equal(actual[["p"]], pf(statistic, 1, df.residual(full), lower.tail = FALSE), tolerance = 1e-9)
      expect_equal(actual[["df2"]], df.residual(full))
    }
  }
})
