if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(220)
  d = expand.grid(g = 1:6, x = -1:1, rep = 1:12)
  d$w = runif(nrow(d), .5, 4)
  d$y = rbinom(nrow(d), 1, -expm1(-exp(-1 + .8 * d$x + d$g / 6)))
  d$cluster = d$rep
  ctl = dbglm_control(glm_max_iter = 1L, tolerate_max_iter = TRUE)
  defer = dbglm_control(glm_max_iter = 1L, tolerate_max_iter = TRUE, compute_inference = FALSE)
  for (formula in list(y ~ x, y ~ x | g)) {
    fit = dbglm(formula, data = d, family = binomial("cloglog"), weights = ~w,
      vcov = ~cluster, control = ctl, retain = "inference")
    fast = dbglm(formula, data = d, family = binomial("cloglog"), weights = ~w,
      vcov = ~cluster, control = defer, retain = "inference")
    expect_false(fast$converged)
    expect_identical(fast$stop_reason, "max_iter")
    expect_equal(coef(fit), coef(fast), tolerance = 1e-12)
    expect_equal(predict(fit, type = "link"), predict(fast, type = "link"), tolerance = 1e-12)
    expect_equal(fitted(fit), fitted(fast), tolerance = 1e-12)
    expect_equal(deviance(fit), deviance(fast))
    expect_true(all(is.na(vcov(fast))))
    expect_true(all(is.finite(vcov(fit))))
    expect_identical(fast$inference_result$reason, "not_requested")
    expect_true(fit$inference_result$conditional_on_iterate)
    expect_true(is.null(fast$retained_inference$H))
    expect_length(fast$retained_inference$design, 0L)
    expect_error(vcov(fast, vcov = "hc1"), "compute_inference=FALSE")
    expect_error(summary(fast, vcov = "hc1"), "compute_inference=FALSE")
    expect_equal(generics::augment(fast)$.fitted, fitted(fast))
    close(fit)
    close(fast)
  }
  expect_error(dbglm_control(compute_inference = NA), "compute_inference")
  expect_error(dbglm(y ~ x | g, data = d, family = binomial("cloglog"),
    control = dbglm_control(glm_max_iter = 1L, compute_inference = FALSE)), "did not converge")
}
