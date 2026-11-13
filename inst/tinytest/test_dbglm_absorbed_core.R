local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = DBI::dbListTables(conn)
  set.seed(11143)
  d = expand.grid(g = 1:5, h = 1:4, replicate = 1:8)
  d[["x"]] = rnorm(nrow(d))
  d[["w"]] = runif(nrow(d), .3, 4)
  d[["o"]] = runif(nrow(d), -.2, .2)
  d[["cl"]] = sample(1:12, nrow(d), TRUE)
  d[["absorbed"]] = d[["g"]] - d[["h"]]
  d[["y"]] = rbinom(nrow(d), 1, plogis(.5 * d[["x"]] + .1 * d[["absorbed"]] + d[["o"]]))
  # Add a singleton. Its high weight must not determine the fitted scale.
  extra = d[1L, ]
  extra[["g"]] = 6L
  extra[["w"]] = 1e8
  data = rbind(d, extra)
  fit = dbglm(y ~ x + absorbed | g + h, data = data, conn = conn,
    family = binomial(), weights = ~w, offset = ~o, vcov = ~cl,
    ssc = "full", start = c(x = .01), null_model = TRUE)
  ref = suppressWarnings(glm(y ~ x + factor(g) + factor(h), d, family = binomial(),
    weights = w, offset = o, control = glm.control(epsilon = 1e-13)))
  X = model.matrix(ref)
  mu = fitted(ref)
  W = d[["w"]] * mu * (1 - mu)
  scores = rowsum(X * (d[["w"]] * (d[["y"]] - mu)), d[["cl"]])
  bread = solve(crossprod(X, W * X))
  N = nrow(d)
  G = nrow(scores)
  V = bread %*% crossprod(scores) %*% bread * G / (G - 1) * (N - 1) / (N - ncol(X))
  expect_identical(fit[["collin.var"]], "absorbed")
  expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-8)
  expect_equal(vcov(fit)["x", "x"], V["x", "x"], tolerance = 1e-8)
  expect_equal(fit[["deviance"]], deviance(ref), tolerance = 1e-10)
  expect_equal(nobs(fit), N)
  expect_equal(fit[["nobs_pruned"]], 1)
  expect_equal(fit[["exclusions"]][["singleton"]], 1)
  expect_true(fit[["converged"]] && fit[["fe_score_max"]] <= 1e-8)
  expect_true(fit[["null_diagnostics"]][["fit"]][["available"]])
  expect_error(predict(fit, d), "retained fixed effects")
  expect_error(dbglm(y ~ x | g + h, data = data, conn = conn,
    family = binomial(), fixef_rm = "none"), "separation remains")
  expect_identical(DBI::dbListTables(conn), initial)

  tuple = dbglm(y ~ x | g^h, data = d, family = gaussian(), weights = ~w,
    offset = ~o, fixef_rm = "none")
  oracle = lm(y ~ x + interaction(g, h), d, weights = w, offset = o)
  expect_equal(coef(tuple)[["x"]], coef(oracle)[["x"]], tolerance = 1e-10)
  expect_equal(unname(tuple[["n_fe_levels"]]), 20)
  only = dbglm(y ~ 1 | g, data = d, family = binomial())
  expect_equal(length(coef(only)), 0)
})
