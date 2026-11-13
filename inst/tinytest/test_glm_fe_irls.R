local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11142)
  d = expand.grid(g = 1:4, h = 1:5, replicate = 1:8)
  n = nrow(d)
  d[["__dbglm_a"]] = runif(n, .2, 2)
  d[["__dbglm_o"]] = runif(n, -.2, .2)
  d[["x"]] = rnorm(n)
  d[["absorbed"]] = cos(d[["g"]]) + sin(d[["h"]])
  eta = .4 * d[["x"]] + .35 * d[["absorbed"]] + d[["__dbglm_o"]]
  fe = c("g", "h")
  sizes = c(4, 5)
  effects = lapply(sizes, numeric)
  control = dbglm_control(score_tol = 1e-9)
  for (family in list(gaussian(), poisson(), binomial("logit"), binomial("cloglog"))) {
    mu = family[["linkinv"]](eta)
    d[["y"]] = switch(family[["family"]], gaussian = rnorm(n, mu),
      poisson = rpois(n, mu), binomial = rbinom(n, 1, mu))
    d[["__dbglm_y"]] = d[["y"]]
    DBI::dbWriteTable(conn, "glm_fe_base", d, overwrite = TRUE)
    dropped = dbreg:::dbglm_fe_irls(conn, "glm_fe_base", c("x", "absorbed"),
      fe, sizes, n, c(0, 0), effects, family, control)
    expect_identical(dropped[["restart"]], 1L)
    fit = dbreg:::dbglm_fe_irls(conn, "glm_fe_base", "x",
      fe, sizes, n, 0, effects, family, control)
    ref = suppressWarnings(glm(y ~ x + factor(g) + factor(h), data = d, family = family,
      weights = d[["__dbglm_a"]], offset = d[["__dbglm_o"]],
      control = glm.control(epsilon = 1e-13, maxit = 100L)))
    expect_true(fit[["converged"]])
    expect_true(max(fit[["score_max"]], fit[["fe_score_max"]]) <= 1e-9)
    expect_equal(fit[["beta"]], unname(coef(ref)[["x"]]), tolerance = 1e-7)
    expect_equal(fit[["deviance"]], deviance(ref), tolerance = 1e-10)
  }
  fit = dbreg:::dbglm_fe_irls(conn, "glm_fe_base", character(), fe, sizes,
    n, numeric(), effects, binomial("cloglog"), control)
  expect_true(fit[["converged"]])
  expect_equal(length(fit[["beta"]]), 0)
})
