# A fitted raw basis must be identified and reconstruct its own predictions.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(51)
  d = data.frame(g = factor(rep(letters[1:3], each = 20)))
  d$x = as.numeric(d$g == "a")
  for (family in list(gaussian(), poisson(), binomial())) {
    d$y = switch(family$family,
      gaussian = 10 + as.integer(d$g) + rnorm(60, sd = .1),
      poisson = rep(c(1, 2, 3, 4), 15), binomial = rep(c(0, 1, 0, 1, 1), 12))
    for (fml in list(y ~ 0 + x + g, y ~ 0 + g + x)) {
      for (strategy in c("moments", "compress")) {
        fit = dbglm(fml, data = d, conn = conn, family = family,
          strategy = strategy, retain = "inference", ssc = list(K.adj = FALSE))
        X = model.matrix(fml, d)[, names(coef(fit)), drop = FALSE]
        expect_equal(qr(X)$rank, ncol(X))
        expected = glm.fit(X, d$y, family = family,
          control = glm.control(epsilon = 1e-12))
        expect_true(fit$converged)
        expect_equal(unname(coef(fit)), unname(coef(expected)), tolerance = 1e-8)
        expect_equal(as.numeric(X %*% coef(fit)), predict(fit, d), tolerance = 1e-10)
        expect_equal(predict(fit, d, type = "response"), fitted(fit), tolerance = 1e-10)
        phi = if (family$family == "gaussian") sum(expected$residuals^2) / (nrow(d) - ncol(X)) else 1
        fisher = family$mu.eta(expected$linear.predictors)^2 / family$variance(expected$fitted.values)
        bread = solve(crossprod(X, X * fisher))
        expect_equal(as.numeric(vcov(fit)), as.numeric(bread * phi), tolerance = 1e-8)
        score = X * (family$mu.eta(expected$linear.predictors) *
          (d$y - expected$fitted.values) / family$variance(expected$fitted.values))
        robust = bread %*% crossprod(score) %*% bread
        expect_equal(as.numeric(vcov(fit, vcov = "hc1")), as.numeric(robust), tolerance = 1e-8)
        close(fit)
      }
    }
  }
})

# Rescaling an absorbed structural predictor must preserve convergence, fitted
# probabilities and inference, including at units where raw Newton steps are big.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(442)
  d = data.frame(g = rep(1:5, each = 30), z = rnorm(150))
  d$y = rbinom(150, 1, plogis(.2 + .8 * d$z + .05 * d$g))
  for (link in c("logit", "probit", "cloglog")) {
    for (method in c("auto", "sql")) {
      control = dbglm_control(projection_method = method)
      base = dbglm(y ~ z | g, data = d, conn = conn,
        family = binomial(link), control = control, retain = "inference")
      for (scale in c(1e-15, 1e-8)) {
        d$x = d$z * scale
        fit = dbglm(y ~ x | g, data = d, conn = conn,
          family = binomial(link), control = control, retain = "inference")
        expect_true(fit$converged)
        expect_equal(unname(coef(fit) * scale), unname(coef(base)), tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit)) * scale^2, as.numeric(vcov(base)), tolerance = 1e-8)
        expect_equal(fitted(fit), fitted(base), tolerance = 1e-8)
        close(fit)
      }
      close(base)
    }
  }
})

# Starting close to the boundary must not turn complete separation into a finite MLE.
local({
  d = data.frame(x = rep(c(-1, 1), each = 20), y = rep(c(0, 1), each = 20))
  expect_error(dbglm(y ~ x, data = d, family = binomial(), start = c(0, 30),
    control = dbglm_control(score_tol = 1e-6)), "probability boundary")
})
