local({
  set.seed(11139)
  n = 200L
  d = data.frame(x = rnorm(n), constant = 1,
    `__dbglm_o` = runif(n, -.25, .25), `__dbglm_a` = runif(n, .2, 1),
    check.names = FALSE)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  controls = list(glm_max_iter = 100L, step_max_halvings = 25L,
    glm_tol = 1e-10, score_tol = 1e-8)
  d[["__dbglm_cluster"]] = rep(1:20, length.out = n)
  ssc = dbreg:::normalize_ssc("full")
  X = model.matrix(~ x, d)
  eta = .2 + .5 * d[["x"]] + d[["__dbglm_o"]]
  for (family in list(gaussian(), poisson(), binomial("logit"),
      binomial("probit"), binomial("cloglog"))) {
    mu = family[["linkinv"]](eta)
    y = switch(family[["family"]], gaussian = rnorm(n, mu),
      poisson = rpois(n, mu), binomial = rbinom(n, 1, mu))
    d[["__dbglm_y"]] = y
    DBI::dbWriteTable(conn, "glm_solver", d, overwrite = TRUE)
    ref = suppressWarnings(stats::glm.fit(X, y, weights = d[["__dbglm_a"]],
      offset = d[["__dbglm_o"]], family = family,
      control = glm.control(epsilon = 1e-13, maxit = 100L)))
    fit = dbreg:::dbglm_irls(conn, "glm_solver", c("constant", "x"),
      c(0, 0), family, controls)
    expect_true(ref[["converged"]])
    expect_equal(fit[["beta"]], unname(ref[["coefficients"]]), tolerance = 1e-7)
    expect_equal(fit[["deviance"]], ref[["deviance"]], tolerance = 1e-10)
    expect_true(fit[["score_max"]] <= controls[["score_tol"]])
    expect_true(all(diff(fit[["deviance_history"]]) < 1e-10))
    eta_fit = as.numeric(X %*% fit[["beta"]]) + d[["__dbglm_o"]]
    mu_fit = family[["linkinv"]](eta_fit)
    dmu = family[["mu.eta"]](eta_fit)
    variance = family[["variance"]](mu_fit)
    bread = solve(crossprod(X, (d[["__dbglm_a"]] * dmu^2 / variance) * X))
    scores = X * (d[["__dbglm_a"]] * dmu * (y - mu_fit) / variance)
    grouped = rowsum(scores, d[["__dbglm_cluster"]])
    p = ncol(X)
    G = nrow(grouped)
    reference = list(iid = bread * (n - 1) / (n - p),
      hc1 = bread %*% crossprod(scores) %*% bread * n / (n - p),
      cluster = bread %*% crossprod(grouped) %*% bread *
        (n - 1) / (n - p) * G / (G - 1))
    for (kind in names(reference)) {
      V = dbreg:::dbglm_covariance(conn, "glm_solver", c("constant", "x"),
        fit, family, kind, ssc, n)[["matrix"]]
      expect_equal(V, unname(reference[[kind]]), tolerance = 1e-10)
    }
  }
  controls[["glm_max_iter"]] = 1L
  controls[["tolerate_max_iter"]] = TRUE
  fit = dbreg:::dbglm_irls(conn, "glm_solver", c("constant", "x"),
    c(0, 0), binomial("cloglog"), controls)
  expect_false(fit[["converged"]])
  expect_identical(fit[["stop_reason"]], "max_iter")
})
