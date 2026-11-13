local({
  prior = options(dbreg.ap_tol = 2e-7)
  on.exit(options(prior))
  expect_equal(dbglm_control()$fe_tol, 2e-7)
  options(prior)
  for (value in list(0, -1, NA_real_, Inf, c(1, 2), "small")) {
    expect_error(dbglm_control(fe_tol = value), "fe_tol")
  }
  for (value in list(0, 1.5, Inf, NA_real_, 2^32)) {
    expect_error(dbglm_control(fe_max_iter = value), "fe_max_iter")
  }
  for (value in list(1, NA, c(TRUE, FALSE), NULL)) {
    expect_error(dbglm_control(tighten_fe_tol = value), "tighten_fe_tol")
  }
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  set.seed(11166)
  d = expand.grid(g = 1:7, h = 1:5, k = 1:3, rep = 1:3)
  d$x = rnorm(nrow(d))
  d$w = runif(nrow(d), .4, 3)
  d$o = sin(d$k)/10
  d$y = rpois(nrow(d), exp(.2 + .3*d$x + cos(d$g)/4 + sin(d$h)/5 + d$o))
  dense = glm(y ~ x + factor(g) + factor(h), data = d, family = poisson(),
    weights = w, offset = o, control = glm.control(epsilon = 1e-14))
  for (method in c("auto", "sql")) {
    control = dbglm_control(projection_method = method, fe_tol = 1e-4,
      score_tol = 1e-9, fe_max_iter = 300L)
    fit = dbglm(y ~ x | g + h, data = d, conn = conn, weights = ~w, offset = ~o,
      family = poisson(), control = control, retain = "prediction")
    expect_equal(fit$projection_tolerance, 1e-11)
    expect_identical(fit$projection_max_iter, 300L)
    expect_true(fit$converged && max(fit$score_max, fit$fe_score_max) <= 1e-9)
    expect_equal(unname(coef(fit)), unname(coef(dense)["x"]), tolerance = 1e-8)
    expect_equal(predict(fit, d, type = "link"), unname(predict(dense, type = "link")), tolerance = 1e-8)
    close(fit)
    control$tighten_fe_tol = FALSE
    control$compute_inference = FALSE
    control$glm_max_iter = 1L
    control$tolerate_max_iter = TRUE
    fit = dbglm(y ~ x | g + h, data = d, conn = conn, weights = ~w, offset = ~o,
      family = poisson(), control = control)
    expect_equal(fit$projection_tolerance, 1e-4)
    expect_false(fit$converged)
    expect_identical(fit$stop_reason, "max_iter")
    expect_true(max(fit$score_max, fit$fe_score_max) > control$score_tol)
    expect_true(all(is.na(vcov(fit))))
    close(fit)
    expect_identical(DBI::dbListTables(conn), character())
  }
  expect_error(dbglm(y ~ x | g + h + k, data = d, conn = conn, family = poisson(),
    weights = ~w, offset = ~o, control = dbglm_control(projection_method = "sql", fe_max_iter = 1L)),
    "iteration cap|broke down|stalled")
  expect_identical(DBI::dbListTables(conn), character())
})
