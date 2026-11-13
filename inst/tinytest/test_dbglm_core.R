local({
  set.seed(11140)
  n = 200L
  d = data.frame(x = rnorm(n), g = factor(rep(letters[1:3], length.out = n)),
    w = runif(n, .1, 3), offset = rnorm(n, sd = .1), cl = rep(1:10, length.out = n))
  d[["z"]] = 2 * d[["x"]]
  d[["use"]] = seq_len(n) %% 7L != 0L
  d[["w"]][1L] = 0
  d[["x"]][2L] = NA_real_
  X = model.matrix(~ x + g, transform(d, x = ifelse(is.na(x), 0, x)))
  eta = as.numeric(X %*% c(.2, .3, .4, -.2)) + d[["offset"]]
  keep = d[["use"]] & d[["w"]] > 0 & !is.na(d[["x"]])
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  for (family in list(gaussian(), poisson(), binomial("logit"),
      binomial("probit"), binomial("cloglog"))) {
    mu = family[["linkinv"]](eta)
    d[["y"]] = switch(family[["family"]], gaussian = rnorm(n, mu),
      poisson = rpois(n, mu), binomial = rbinom(n, 1, mu))
    ref = suppressWarnings(glm(y ~ x + g, data = d[keep, ], family = family,
      weights = w, offset = offset, control = glm.control(epsilon = 1e-15, maxit = 1000L)))
    fit = dbglm(y ~ x + z + g, data = d, conn = conn, family = family,
      weights = ~w, offset = ~offset, subset = ~use, vcov = ~cl, null_model = TRUE)
    expect_identical(names(coef(fit)), names(coef(ref)))
    expect_identical(fit[["collin.var"]], "z")
    expect_true(max(abs(coef(fit) - coef(ref))) < 1e-7)
    expect_equal(fit[["deviance"]], deviance(ref), tolerance = 1e-10)
    likelihood = switch(family[["family"]],
      gaussian = as.numeric(logLik(ref)) - 1,
      poisson = as.numeric(logLik(ref)),
      binomial = sum(round(d[["w"]][keep]) *
        dbinom(d[["y"]][keep], 1, predict(fit, d[keep, ], type = "response"), log = TRUE)))
    expect_equal(as.numeric(logLik(fit)), likelihood, tolerance = 1e-10)
    null_ref = suppressWarnings(glm(y ~ 1, data = d[keep, ], family = family,
      weights = w, offset = offset, control = glm.control(epsilon = 1e-15, maxit = 1000L)))
    null = fit[["null_diagnostics"]]
    expect_true(null[["fit"]][["available"]])
    expect_equal(null[["nobs"]], sum(keep))
    expect_equal(null[["reported"]][["logLik"]],
      as.numeric(logLik(null_ref)) - as.numeric(family[["family"]] == "gaussian"),
      tolerance = 1e-7)
    expect_equal(null[["reported"]][["ratio"]],
      1 - as.numeric(logLik(fit)) / null[["reported"]][["logLik"]])
    expect_equal(nobs(fit), sum(keep))
    expect_true(fit[["converged"]] && fit[["score_max"]] <= 1e-8)
    expect_equal(predict(fit, d[keep, ], type = "response"), unname(predict(ref, type = "response")),
      tolerance = 1e-7)
    expect_true(max(abs(predict(fit, d[keep, ]) - unname(predict(ref, type = "link")))) < 1e-7)
    expect_true(is.na(predict(fit, d[2L, ])))
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  d[["y"]] = rnorm(n)
  fit = dbglm(y ~ 0 + x, data = d, start = c(x = .01), null_model = TRUE)
  expect_identical(fit[["null_diagnostics"]][["fit"]][["method"]], "fixed_offset")
  expect_equal(coef(fit), coef(lm(y ~ 0 + x, d)), tolerance = 1e-10)
  expect_error(dbglm(y ~ x, data = d, start = c(a = 0, b = 0)), "Named starts")
  expect_error(dbglm(y ~ x | g | z, data = d), "structural formula")
  expect_error(dbglm(y ~ x, data = d, family = poisson()), "must be nonnegative")
  previous = options(dbreg.projection_max_entries = 1L)
  expect_error(dbglm(y ~ x, data = d, conn = conn), "information exceeds")
  options(previous)
  d[["w"]][3L] = Inf
  expect_error(dbglm(y ~ x, data = d, weights = ~w, conn = conn), "finite")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
