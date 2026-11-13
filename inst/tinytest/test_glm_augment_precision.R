# Augmentation must use the response residual saved before restoring a large
# Gaussian location, just as residuals() does.
local({
  set.seed(5)
  d = data.frame(x = rnorm(200), g = rep(1:10, each = 20), o = rep(c(-.1, .1), 100))
  d[["y"]] = 1e12 + .002*d[["x"]] + rnorm(200, sd = .0005)
  for (fml in list(y ~ x, y ~ x | g)) for (strategy in c("moments", "compress")) {
    for (with_offset in c(FALSE, TRUE)) {
      data = d
      data[["target"]] = data[["y"]] - 1e12
      if (with_offset) data[["y"]] = data[["y"]] + data[["o"]]
      # Center the actually stored response/offset difference, including its
      # representable input rounding, rather than the simulation's unrounded y.
      data[["target"]] = (data[["y"]] - if (with_offset) data[["o"]] else 0) - 1e12
      ref = lm(if (length(Formula::Formula(fml))[[2L]] > 1L) target ~ x + factor(g) else target ~ x,
        data = data)
      fit = dbglm(fml, data = data, strategy = strategy, retain = "inference",
        offset = if (with_offset) "o" else NULL)
      actual = generics::augment(fit)
      expect_equal(actual[[".resid"]], residuals(fit), tolerance = 1e-14)
      expect_equal(actual[[".resid"]], unname(residuals(ref)), tolerance = 1e-10)
      expect_equal(generics::augment(fit, data = data)[[".resid"]], actual[[".resid"]], tolerance = 1e-14)
      close(fit)
    }
  }
})
