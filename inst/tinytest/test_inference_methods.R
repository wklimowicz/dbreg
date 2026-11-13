if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(156)
  n = 160L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = rep(1:8, length.out = n),
    cl = sample(1:11, n, replace = TRUE))
  d$y = d$x - d$z + rnorm(n, sd = exp(d$x / 4))
  d$count = rpois(n, exp(d$x / 3))
  fits = list(dbreg(y ~ x + z | g, data = d, retain = "inference", strategy = "demean"))
  if (exists("dbglm", mode = "function")) {
    fits = c(fits, list(dbglm(count ~ x + z, data = d, family = poisson(), retain = "inference")))
  }
  for (fit in fits) {
    original = vcov(fit)
    V = vcov(fit, cluster = "cl", ssc = list(t.df = 4))
    sm = summary(fit, cluster = "cl", ssc = list(t.df = 4))
    expect_true(inherits(sm, "summary.dbreg"))
    expect_identical(vcov(fit), original)
    expect_equal(vcov(sm), V)
    expect_equal(coef(sm), coef(fit))
    expect_equal(vcov(sm, ssc = list(t.df = 7)), vcov(fit, cluster = "cl", ssc = list(t.df = 7)))
    se = sqrt(Matrix::diag(V))
    normal = inherits(fit, "dbglm")
    q = if (normal) qnorm(.975) else qt(.975, 4)
    ci = confint(fit, cluster = "cl", ssc = list(t.df = 4))
    expect_equal(ci[, 2] - ci[, 1], 2 * q * se, tolerance = 1e-12)
    expect_equal(confint(sm), ci)
    tidied = tidy(fit, conf.int = TRUE, cluster = "cl", ssc = list(t.df = 4))
    expect_equal(tidied$std.error, unname(se))
    expect_equal(tidied$conf.low, unname(ci[, 1]))
    t = coef(fit) / se
    expect_equal(tidied$p.value, unname(2 * if (normal) pnorm(-abs(t)) else pt(-abs(t), 4)))
    if (requireNamespace("fixest", quietly = TRUE)) {
      table = fixest::coeftable(sm, keep = c("^x$", "^z$"), order = "^z$")
      expect_identical(rownames(table), c("z", "x"))
      expect_identical(colnames(table)[3L], if (normal) "z value" else "t value")
      expect_equal(fixest::se(sm), se)
      expect_equal(fixest::tstat(sm), t)
      expect_equal(fixest::pvalue(sm), setNames(tidied$p.value, tidied$term))
    }
    printed = capture.output(print(sm))
    expect_true(any(grepl(if (normal) "GLM estimation" else "OLS estimation", printed)))
    close(fit)
    expect_equal(vcov(sm), V)
  }
}
