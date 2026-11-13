if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(11161)
  d = data.frame(x = rnorm(120), z = rnorm(120), g = rep(1:6, each = 20),
    h = rep(1:10, 12))
  d$copy = d$x
  d$y = .5 * d$x + .4 * d$z + sin(d$g) + rnorm(nrow(d))
  d$count = rpois(nrow(d), exp(.2 + .3 * d$x))
  fit = dbreg(y ~ x + z + copy | g, data = d, strategy = "demean", retain = "all")
  expect_identical(names(coef(fit, collin = TRUE)), c("x", "z", "copy"))
  expect_true(is.na(coef(fit, collin = TRUE)[["copy"]]))
  V = vcov(fit, vcov = "hc1")
  df = attr(V, "df.t")
  supplied = as.matrix(V)[2:1, 2:1]
  attr(supplied, "df.t") = df
  coefficients = tidy(fit, vcov = supplied, conf.int = TRUE, collin = TRUE)
  expect_identical(coefficients$term, c("x", "z", "copy"))
  expect_true(all(is.na(coefficients[3L, -1L])))
  expect_equal(coefficients$std.error[1:2], unname(sqrt(diag(V))))
  expect_equal(coefficients$conf.low[1:2], unname(confint(fit, vcov = "hc1")[, 1L]))
  expect_error(tidy(fit, vcov = matrix(1, 1, 1)), "symmetric covariance")
  expect_error(tidy(fit, vcov = supplied, cluster = ~g), "cannot be combined")
  two = suppressWarnings(summary(fit, vcov = ~g + h))
  expect_equal(nrow(glance(two)), 1L)
  expect_identical(glance(two)$vcov.type, "twoway")

  if (requireNamespace("insight", quietly = TRUE)) {
    expect_equal(insight::get_varcov(fit, vcov = "HC1"), V)
    expect_equal(insight::get_varcov(fit, vcov = "vcovCL", vcov_args = list(cluster = ~g)),
      vcov(fit, cluster = ~g))
    expect_error(insight::get_varcov(fit, vcov = "HC3"), "Unsupported covariance")
  }
  if (requireNamespace("modelsummary", quietly = TRUE) && requireNamespace("broom", quietly = TRUE)) {
    for (request in list("HC1", ~g)) {
      estimates = modelsummary::get_estimates(fit, vcov = request)
      expect_equal(estimates$std.error, unname(sqrt(diag(vcov(fit, vcov = request)))))
      expect_equal(estimates$conf.low, unname(confint(fit, vcov = request)[, 1L]))
    }
  }
  single = dbreg(y ~ x | g, data = d, strategy = "demean")
  if (requireNamespace("fixest", quietly = TRUE)) {
    table = fixest::coeftable(single)
    expect_identical(names(fixest::se(single)), "x")
    expect_equal(fixest::se(single), sqrt(diag(vcov(single))))
    expect_equal(fixest::tstat(single), setNames(as.numeric(table[, 3L]), "x"))
    expect_equal(fixest::pvalue(single), setNames(as.numeric(table[, 4L]), "x"))
  }
  close(single)
  close(fit)
  fit = dbglm(count ~ x, data = d, family = poisson(), null_model = TRUE)
  dense = glm(count ~ x, data = d, family = poisson())
  dense_null = glm(count ~ 1, data = d, family = poisson())
  summary = glance(fit)
  expect_equal(nrow(summary), 1L)
  expect_equal(summary$logLik, as.numeric(logLik(dense)), tolerance = 1e-10)
  expect_equal(summary$pseudo.r.squared, 1 - as.numeric(logLik(dense) / logLik(dense_null)), tolerance = 1e-10)
  expect_equal(unname(gof(fit)[["null.logLik"]]), as.numeric(logLik(dense_null)), tolerance = 1e-10)
}
