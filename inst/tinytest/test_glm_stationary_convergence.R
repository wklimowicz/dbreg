# Controlled objective roundoff proves driver stopping behavior; these fixtures
# are not substitutes for a native distributed convergence certificate.
local({
  ns = asNamespace("dbreg")
  probe = function(absorbed, previous = 1e-12, current = 1e-12,
                   previous_fe = 1e-12, current_fe = 1e-12,
                   boundary = FALSE, invalid = FALSE) {
    driver = get(if (absorbed) "dbglm_fe_irls" else "dbglm_irls", ns)
    e = new.env(parent = environment(driver))
    calls = 0L
    state = function(...) {
      calls <<- calls + 1L
      starting = calls == 1L
      score = if (starting) previous else current
      list(bad = if (!starting && invalid) 1L else 0L,
        mu_min = .5, mu_max = if (!starting && boundary) Inf else .5,
        deviance = if (!starting && calls %% 2L == 0L) 10.01 else 10,
        H = matrix(1, 1, 1), score = score, score_max = abs(score),
        working = as.character(calls))
    }
    e$dbglm_fe_trial = state
    e$dbglm_evaluate = state
    e$dbglm_fe_step = function(...) list(H = matrix(1, 1, 1), score = 0,
      projection = list(matrix(0, 1, 2)), residual_source = "mock")
    e$dbglm_effect_delta = function(...) list(0)
    e$dbglm_fe_group_score = function(conn, working, ...) {
      if (working == "1") previous_fe else current_fe
    }
    environment(driver) = e
    control = dbglm_control(glm_max_iter = 1L, step_max_halvings = 1L,
      tolerate_max_iter = TRUE)
    tryCatch(if (absorbed) driver(NULL, "mock", "x", "g", 1L, 10L, 0,
      list(0), list(link = "identity"), control) else
        driver(NULL, "mock", "x", 0, list(link = "identity"), control),
      error = identity)
  }
  for (absorbed in c(FALSE, TRUE)) {
    fit = probe(absorbed)
    expect_true(isTRUE(fit$converged))
    expect_equal(fit$iterations, 1L)
    expect_equal(fit$step_halvings, 1L)
    expect_equal(fit$deviance_history, c(10, 10))
    for (previous in c(1e-3, Inf)) {
      fit = probe(absorbed, previous = previous)
      expect_false(isTRUE(fit$converged))
    }
    for (current in c(1e-3, Inf, NA_real_)) {
      fit = probe(absorbed, current = current)
      expect_false(isTRUE(fit$converged))
    }
    expect_false(isTRUE(probe(absorbed, boundary = TRUE)$converged))
    expect_false(isTRUE(probe(absorbed, invalid = TRUE)$converged))
  }
  expect_false(isTRUE(probe(TRUE, previous_fe = 1e-3)$converged))
  expect_false(isTRUE(probe(TRUE, current_fe = 1e-3)$converged))
})

local({
  converged = dbreg:::dbglm_converged_step
  norm = dbreg:::dbglm_score_norm
  control = dbglm_control()
  expect_true(converged(.5, 0, "logit", 1e-12, 1e-12, control))
  expect_true(converged(1, 0, "logit", 1e-12, Inf, control))
  expect_false(converged(.5, 0, "logit", 1e-12, 1e-3, control))
  expect_false(converged(.5, 1, "logit", 1e-12, 1e-12, control))
  expect_false(converged(.5, NA_real_, "logit", 1e-12, 1e-12, control))
  for (bad in c(NA_real_, Inf, -1, 1e-3)) {
    expect_false(converged(.5, 0, "identity", bad, 1e-12, control))
    expect_false(converged(.5, 0, "identity", 1e-12, bad, control))
  }
  expect_equal(norm(c(1, 2), diag(c(4, 16))), .5)
  expect_equal(norm(numeric(), matrix(numeric(), 0, 0)), 0)
  for (H in list(matrix(Inf, 1, 1), matrix(0, 1, 1),
      matrix(-1, 1, 1), matrix(1, 2, 2), diag(c(1, NaN)))) {
    expect_identical(norm(rep(0, nrow(H)), H), Inf)
  }
  expect_identical(norm(0, diag(2)), Inf)
  expect_identical(norm(NA_real_, matrix(1, 1, 1)), Inf)
})
