# Promotion must happen inside interactions/powers, before their first product.
local({
  set.seed(421)
  d = data.frame(x = 50000L+seq_len(200), z = 60000L+sample.int(200), y = rnorm(200))
  for (formula in list(y ~ x:z, y ~ I(x^2))) {
    reference = stats::lm(formula, d)
    for (strategy in c("moments", "compress")) {
      fit = dbreg(formula, data = d, strategy = strategy)
      expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-7)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-7)
    }
    gaussian = dbglm(formula, data = d, family = stats::gaussian())
    expect_equal(unname(coef(gaussian)), unname(coef(reference)), tolerance = 1e-6)
  }
})
