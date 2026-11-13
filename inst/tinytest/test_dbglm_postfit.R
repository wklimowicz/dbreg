if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(917)
  d = expand.grid(g = factor(1:5), x = -1:1, trend = -1:1, rep = 1:12)
  d$w = sample(c(0.5, 2, 8), nrow(d), replace = TRUE)
  d$o = d$trend / 10
  d$y = rpois(nrow(d), exp(0.2 + d$x / 4 + as.integer(d$g) * d$trend / 10 + d$o))
  d[["new cluster"]] = (d$rep + as.integer(d$g)) %% 11
  d$batch = d$rep %% 2
  d$bad = d[["new cluster"]]
  d$bad[1] = NA
  correction = list(K.adj = FALSE, G.adj = FALSE)
  bare = function(V) matrix(as.numeric(V), nrow(V), ncol(V))
  covariance_reference = function(model, columns, groups = NULL) {
    X = model.matrix(model)
    mu = fitted(model)
    bread = solve(crossprod(X, X * (model$prior.weights * mu)))
    scores = X * (model$prior.weights * (model$y - mu))
    if (!is.null(groups)) scores = rowsum(scores, groups)
    V = bread %*% crossprod(scores) %*% bread
    j = match(columns, colnames(X))
    bare(V[j, j, drop = FALSE])
  }
  cases = list(
    list(y ~ x, y ~ x, "compress", "auto"),
    list(y ~ x | g, y ~ x + g, "moments", "auto"),
    list(y ~ x | g[trend], y ~ x + g / trend, "compress", "sql")
  )
  for (case in cases) {
    ref = glm(case[[2]], data = d, family = poisson(), weights = w, offset = o)
    fit = dbglm(case[[1]], data = d, family = poisson(), weights = ~w, offset = ~o,
      retain = "inference", strategy = case[[3]],
      control = dbglm_control(projection_method = case[[4]]))
    columns = names(coef(fit))
    expect_equal(bare(vcov(fit, vcov = "hc1", ssc = correction)),
      covariance_reference(ref, columns), tolerance = 1e-8)
    V = vcov(fit, vcov = ~`new cluster`^batch, ssc = correction)
    expect_equal(bare(V), covariance_reference(ref, columns, interaction(d[["new cluster"]], d$batch)),
      tolerance = 1e-8)
    expect_equal(attr(V, "n_clusters"), length(unique(interaction(d[["new cluster"]], d$batch))))
    expect_equal(vcov(fit, vcov = "iid"), vcov(fit), tolerance = 1e-10)
    expect_error(vcov(fit, cluster = "bad"), "missing or nonfinite")
    expect_error(vcov(fit, cluster = "absent"), "absent from the retained source")
    expect_error(vcov(fit, cluster = c("g", "batch")), "one cluster term")
    close(fit)
    expect_error(vcov(fit, vcov = "hc1"), "no live database session")
  }

  s = list(K.fixef = "nonnested")
  pre = dbglm(y ~ x | g, data = d, family = poisson(), weights = ~w, offset = ~o,
    vcov = ~g, retain = "inference", ssc = s)
  expect_equal(vcov(pre, cluster = "g"), vcov(pre), tolerance = 1e-10)
  expect_equal(vcov(pre, ssc = s), vcov(pre), tolerance = 1e-10)
  close(pre)

  d$binary = rbinom(nrow(d), 1, pnorm(d$x / 3 + d$o))
  for (family in list(gaussian(), binomial("probit"))) {
    formula = if (family$family == "gaussian") y ~ x else binary ~ x
    ref = suppressWarnings(glm(formula, data = d, family = family, weights = w, offset = o))
    fit = dbglm(formula, data = d, family = family, weights = ~w, offset = ~o,
      retain = "inference", strategy = "compress")
    X = model.matrix(ref)
    mu = fitted(ref)
    derivative = family$mu.eta(predict(ref))
    variance = family$variance(mu)
    bread = solve(crossprod(X, X * (d$w * derivative^2 / variance)))
    score = X * (d$w * derivative * (ref$y - mu) / variance)
    expect_equal(bare(vcov(fit, vcov = "hc1", ssc = correction)),
      bare(bread %*% crossprod(score) %*% bread), tolerance = 1e-7)
    expect_equal(vcov(fit, vcov = "iid"), vcov(fit), tolerance = 1e-10)
    close(fit)
  }

  none = dbglm(y ~ x, data = d, family = poisson(), retain = "prediction")
  expect_error(vcov(none, vcov = "hc1"), "retain")
  close(none)
}
