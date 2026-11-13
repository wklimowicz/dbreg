local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  set.seed(318)
  d = data.frame(x = rnorm(240), g = rep(1:12, each = 20))
  d$o = rep(c(-.2, .2), 120)
  d$y = rbinom(nrow(d), 1, plogis(.1 + .3 * d$x + d$o))
  oracle = glm(y ~ x + factor(g), data = d, family = binomial(), offset = o,
    control = glm.control(epsilon = 1e-14, maxit = 1000))
  for (projection in c("auto", "sql")) for (start in c(10, 100, 1000)) {
    fit = tryCatch(dbglm(y ~ x | g, data = d, conn = conn, family = binomial(),
      offset = ~o, start = c(x = start), retain = "all",
      control = dbglm_control(projection_method = projection)), error = identity)
    expect_false(inherits(fit, "error"))
    if (!inherits(fit, "error")) {
      expect_true(fit$converged)
      expect_equal(coef(fit)[["x"]], coef(oracle)[["x"]], tolerance = 1e-7)
      expect_equal(predict(fit, d), unname(predict(oracle)), tolerance = 1e-7)
      expect_true(fit$score_max <= 1e-8 && fit$fe_score_max <= 1e-8)
      expect_true(fit$step_halvings > 0)
      close(fit)
    }
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  set.seed(7312)
  d = expand.grid(g = 1:4, h = 1:3, x = -1:1, trend = -1:1, rep = 1:4)
  d$o = d$h / 20
  d$w = rep(c(.5, 1, 2), length.out = nrow(d))
  eta = .2 + .2 * d$x + .03 * d$g * d$trend + d$o
  cases = list(
    list(binomial("probit"), y ~ x | g + h,
      y ~ x + factor(g) + factor(h), "moments", "auto", "fisher", FALSE),
    list(binomial("cloglog"), y ~ x | g[trend] + h,
      y ~ x + factor(g) / trend + factor(h), "compress", "sql", "newton", FALSE),
    list(binomial("cloglog"), y ~ x | g[trend] + h,
      y ~ x + factor(g) / trend + factor(h), "compress", "sql", "fisher", FALSE),
    list(binomial(), y ~ x | g + h,
      y ~ x + factor(g) + factor(h), "moments", "sql", "fisher", TRUE)
  )
  for (case in cases) {
    family = case[[1L]]
    d$y = rbinom(nrow(d), 1, family$linkinv(eta))
    d$warm = d$o + 1000 * d$x + 100 * d$g
    oracle = suppressWarnings(glm(case[[3L]], data = d, weights = w, offset = o,
      family = family, control = glm.control(epsilon = 1e-14, maxit = 1000)))
    fit = dbglm(case[[2L]], data = d, conn = conn, weights = ~w, offset = ~o,
      family = family, retain = "all", strategy = case[[4L]], start = c(x = 1000),
      etastart = if (case[[7L]]) ~warm else NULL,
      control = dbglm_control(projection_method = case[[5L]], irls = case[[6L]]))
    expect_true(fit$converged && fit$finite_mle)
    expect_equal(coef(fit)[["x"]], coef(oracle)[["x"]], tolerance = 1e-7)
    expect_equal(predict(fit, d), unname(predict(oracle)), tolerance = 1e-7)
    mu = family$linkinv(predict(fit, d))
    derivative = family$mu.eta(predict(fit, d))
    X = model.matrix(oracle)
    info = crossprod(X, X * (d$w * derivative^2 / family$variance(mu)))
    score = crossprod(X, d$w * derivative * (d$y - mu) / family$variance(mu))
    expect_true(max(abs(score) / sqrt(diag(info))) < 1e-7)
    expect_true(fit$step_halvings > 0)
    close(fit)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # The valid default intercept cancels a large known offset. Shrinking toward
  # a zero predictor would destroy that cancellation during recovery.
  d$o = 200
  d$y = rbinom(nrow(d), 1, plogis(.2 + .2 * d$x))
  oracle = glm(y ~ x, data = d, offset = o, family = binomial(),
    control = glm.control(epsilon = 1e-14))
  f = dbglm(y ~ x, data = d, conn = conn, offset = ~o, family = binomial(),
    start = c("(Intercept)" = 1000, x = 1000))
  expect_equal(coef(f), coef(oracle), tolerance = 1e-7)
  expect_true(f$converged && f$step_halvings > 0)
  close(f)
  for (projection in c("auto", "sql")) {
    expect_error(dbglm(y ~ x | g, data = d, conn = conn, family = binomial(),
      offset = ~o, start = c(x = 1e300), retain = "all",
      control = dbglm_control(step_max_halvings = 1L, tolerate_max_iter = TRUE,
        projection_method = projection)), "Initial GLM predictor")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  expect_error(dbglm(y ~ x | g, data = d, conn = conn, family = binomial(),
    start = c(x = NA_real_)), "finite")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})

local({
  # A valid mean can still give an ill-conditioned first projection. Retry only
  # the classified numerical failures; never absorb DBI/capacity/iteration errors.
  ns = asNamespace("dbreg")
  probe = function(reason) {
    driver = get("dbglm_fe_irls", ns)
    e = new.env(parent = environment(driver))
    attempts = 0L
    e$dbglm_fe_trial = function(...) list(bad = 0L, mu_min = .5, mu_max = .5,
      deviance = 10, H = matrix(1, 1, 1), score = 0, working = "mock")
    e$dbglm_fe_step = function(...) {
      attempts <<- attempts + 1L
      if (attempts == 1L) {
        if (reason == "driver") stop("unclassified driver failure")
        dbreg:::projection_abort(reason, "initial projection failure")
      }
      list(H = matrix(1, 1, 1), score = 0,
        projection = list(matrix(0, 1, 2)), residual_source = "mock")
    }
    e$dbglm_effect_delta = function(...) list(0)
    e$dbglm_fe_group_score = function(...) 0
    environment(driver) = e
    value = tryCatch(driver(NULL, "mock", "x", "g", 1L, 10L, 1000, list(0),
      list(link = "identity"), dbglm_control()), error = identity)
    list(value = value, attempts = attempts)
  }
  for (reason in c("certificate_failed", "singular")) {
    r = probe(reason)
    expect_false(inherits(r$value, "error"))
    if (!inherits(r$value, "error")) {
      expect_true(r$value$converged)
      expect_equal(r$value$step_halvings, 1L)
    }
  }
  for (reason in c("driver", "system_limit", "nonconverged")) {
    r = probe(reason)
    expect_true(inherits(r$value, "error"))
    expect_equal(r$attempts, 1L)
  }
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(g = rep(1:2, each = 5), w = 1, score = (1:10) * 1e-6, x = 1:10)
  DBI::dbWriteTable(conn, "working", d)
  resources = dbreg:::fit_resources()
  resources$conn = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  borrowed = manager$create("SELECT 1 AS keep")
  initial = sort(DBI::dbListTables(conn))
  owned = resources$tables
  projection = dbreg:::dbglm_fe_sql_step
  e = new.env(parent = environment(projection))
  e$fe_projection_certificate = function(...) list(max_abs = Inf, max_relative = Inf)
  environment(projection) = e
  expect_error(projection(conn, "working", "x", "g", nrow(d), 1e-10,
    manager, FALSE, FALSE, dbreg:::fe_cg_control()), "certificate")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  expect_identical(resources$tables, owned)
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 1)
})
