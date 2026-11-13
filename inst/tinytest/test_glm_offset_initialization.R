# The default and recovery anchor must remain usable for heterogeneous offsets.
# Independently fit from a known valid start; there is overlap at offset zero.
local({
  d = data.frame(y = c(rep(0, 20), rep(c(0, 1), 10)),
    o = rep(c(-100, 0), each = 20), g = rep(rep(1:2, each = 5), 4))
  for (family in list(binomial(), binomial("probit"), binomial("cloglog"))) {
    for (absorbed in c(FALSE, TRUE)) {
      fml = if (absorbed) y ~ 1 | g else y ~ 1
      ref = suppressWarnings(glm(if (absorbed) y ~ factor(g) else y ~ 1,
        data = d, offset = o, family = family,
        start = if (absorbed) c(0, 0) else 0,
        control = glm.control(epsilon = 1e-14, maxit = 1000L)))
      for (strategy in c("moments", "compress")) {
        fit = dbglm(fml, data = d, offset = "o", family = family, strategy = strategy,
          retain = "all")
        expect_true(fit[["converged"]] && fit[["finite_mle"]])
        expect_equal(fit[["deviance"]], deviance(ref), tolerance = 1e-8)
        expect_equal(predict(fit, d, type = "response"), unname(fitted(ref)), tolerance = 1e-8)
        close(fit)
      }
    }
  }
})
