# Reuse only matching Fisher moments, while comparing actual projected SQL
# scores and information against a separately residualized weighted QR design.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = expand.grid(g = 1:5, h = 1:4, rep = 1:12)
  i = seq_len(nrow(d))
  d$x = sin(i / 7) + d$g / 20
  d$z = cos(i / 11) - d$h / 10
  d$w = rep(c(.5, 1.5, 4), length.out = nrow(d))
  d$score = d$w * (sin(i / 3) + .3 * d$g)
  DBI::dbWriteTable(conn, "working", d)
  Y = cbind(d$score / d$w, d$x, d$z)
  D = model.matrix(~ factor(g) + factor(h), d)
  expected = Y - D %*% qr.solve(D * sqrt(d$w), Y * sqrt(d$w))
  gram = crossprod(expected, expected * d$w)
  state = list(H = crossprod(Y[, -1L], Y[, -1L] * d$w), score = as.numeric(crossprod(Y[, -1L], d$score)))
  state$working = "working"
  state$working_columns = c("x", "z")
  state$information_type = "fisher"
  project = dbreg:::dbglm_fe_step
  e = new.env(parent = environment(project))
  counts = new.env(parent = emptyenv())
  counts$raw = counts$projected = 0L
  e$fe_information_query = function(conn, source, columns, budget, limits, bounds = NULL) {
    if (grepl("glm_within_coordinates", source, fixed = TRUE)) counts$projected = counts$projected + 1L else counts$raw = counts$raw + 1L
    dbreg:::fe_information_query(conn, source, columns, budget, limits)
  }
  e$fe_fwl_information = dbreg:::fe_fwl_information
  environment(e$fe_fwl_information) = e
  e$dbglm_fe_information = dbreg:::dbglm_fe_information
  environment(e$dbglm_fe_information) = e
  environment(project) = e
  for (cached in c(FALSE, TRUE)) {
    counts$raw = counts$projected = 0L
    fit = project(conn, "working", c("x", "z"), c("g", "h"), c(5L, 4L), nrow(d), 1e-10,
      raw_state = if (cached) state else NULL)
    expect_equal(unname(fit$H), unname(gram[-1L, -1L]), tolerance = 1e-9)
    expect_equal(fit$score, as.numeric(gram[-1L, 1L]), tolerance = 1e-9)
    expect_false(fit$information_recomputed)
    expect_identical(fit$reused_working_moments, cached)
    expect_equal(counts$raw, as.integer(!cached))
    expect_equal(counts$projected, 0L)
    expect_true(max(fit$certificate$max_abs, fit$certificate$max_relative) <= 1e-10)
  }
  expect_error(project(conn, "working", c("x", "z"), c("g", "h"), c(5L, 4L), nrow(d), 1e-10,
    raw_state = list(H = matrix(1), score = 1)), "matching GLM working moments")
  for (field in c("working", "working_columns", "information_type")) {
    invalid = state
    invalid[[field]] = switch(field, working = "other_working_state",
      working_columns = c("z", "x"), information_type = "observed")
    expect_error(project(conn, "working", c("x", "z"), c("g", "h"), c(5L, 4L), nrow(d), 1e-10,
      raw_state = invalid), "matching GLM working moments")
  }
  # Nearly absorbed axes require actual-row products, with the same final
  # projection certificate; never use cancellation to decide their rank.
  d$x = 1e4 * d$g + sin(i / 7)
  DBI::dbWriteTable(conn, "working", d, overwrite = TRUE)
  Y = cbind(d$score / d$w, d$x, d$z)
  expected = Y - D %*% qr.solve(D * sqrt(d$w), Y * sqrt(d$w))
  gram = crossprod(expected, expected * d$w)
  fit = project(conn, "working", c("x", "z"), c("g", "h"), c(5L, 4L), nrow(d), 1e-8)
  expect_true(fit$information_recomputed)
  expect_equal(unname(fit$H), unname(gram[-1L, -1L]), tolerance = 1e-7)
  expect_equal(fit$score, as.numeric(gram[-1L, 1L]), tolerance = 1e-7)
  DBI::dbRemoveTable(conn, "working")
  expect_identical(DBI::dbListTables(conn), character())
})

# A plausible diagonal cannot certify an uncertain structural rank. Recompute
# the real Gram before basis selection, and preserve genuinely dependent axes.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(w = c(.5, 1.5, 2, 1), v1 = c(.3, -.1, .2, -.4),
    v2 = c(-1, -1, 1, 1), v3 = c(-1, 1, -1, 1))
  DBI::dbWriteTable(conn, "within_rows", d)
  budget = function() list2env(list(rows = 0L, cells = 0L, upload_bytes = 0L), parent = emptyenv())
  info = function(H, raw = diag(3)) list(cross = raw, within = rbind(c(1, 0, 0), cbind(0, H)))
  calculate = dbreg:::dbglm_fe_information
  for (H in list(matrix(1, 2, 2), matrix(c(1, 1.1, 1.1, 1), 2),
      matrix(c(1, 1 - 1e-10, 1 - 1e-10, 1), 2), matrix(NA_real_, 2, 2))) {
    value = calculate(conn, "within_rows", c("v1", "v2", "v3"), info(H), budget(), dbreg:::projection_limits(), 4L)
    expect_true(value$recomputed)
    expected = crossprod(as.matrix(d[c("v2", "v3")]), as.matrix(d[c("v2", "v3")]) * d$w)
    expect_equal(unname(value$H), unname(expected), tolerance = 1e-12)
    expect_equal(value$score, as.numeric(crossprod(as.matrix(d[c("v2", "v3")]), d$w * d$v1)), tolerance = 1e-12)
    expect_equal(length(dbreg:::dbglm_identified_slopes(value$H, expected, 1e-9)), 2L)
  }
  d$v3 = 2 * d$v2
  DBI::dbWriteTable(conn, "within_rows", d, overwrite = TRUE)
  value = calculate(conn, "within_rows", c("v1", "v2", "v3"), info(matrix(1, 2, 2)), budget(), dbreg:::projection_limits(), 4L)
  expect_true(value$recomputed)
  expect_equal(length(dbreg:::dbglm_identified_slopes(value$H, value$H, 1e-9)), 1L)
  # Count, transfer and pure-FE cases use the actual rows even with cached raw
  # moments; no approximate moment state may manufacture a certificate.
  expect_error(calculate(conn, "within_rows", c("v1", "v2", "v3"), info(diag(2)), budget(), dbreg:::projection_limits(), 5L), "score certificate")
  limits = dbreg:::projection_limits()
  limits$transfer_cells = 2L
  expect_error(calculate(conn, "absent", c("v1", "v2", "v3"), info(diag(2)), budget(), limits, 4L), "scores exceed")
  pure = calculate(conn, "within_rows", "v1", list(cross = matrix(1), within = matrix(1)), budget(), dbreg:::projection_limits(), 4L)
  expect_equal(dim(pure$H), c(0L, 0L))
  expect_identical(pure$score, numeric())
  expect_false(pure$recomputed)
  DBI::dbRemoveTable(conn, "within_rows")
  expect_identical(DBI::dbListTables(conn), character())
})

# Fisher raw moments describe a different weight system from a noncanonical
# observed-Newton direction. Only Fisher projections receive the matching cache.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(8821)
  d = expand.grid(g = 1:4, x = c(-1, -.3, .4, 1), rep = 1:15)
  d$y = rbinom(nrow(d), 1, binomial("cloglog")$linkinv(-.6 + .25 * d$x + d$g / 10))
  oracle = glm(y ~ x + factor(g), data = d, family = binomial("cloglog"),
    control = glm.control(epsilon = 1e-14, maxit = 1000))
  DBI::dbWriteTable(conn, "base", data.frame("__dbglm_y" = d$y, "__dbglm_a" = 1,
    "__dbglm_o" = 0, x = d$x, g = d$g, check.names = FALSE))
  driver = dbreg:::dbglm_fe_irls
  e = new.env(parent = environment(driver))
  cache = new.env(parent = emptyenv())
  cache$observed = cache$fisher = 0L
  e$dbglm_fe_step = function(conn, working, ..., raw_state = NULL) {
    if (is.null(raw_state)) cache$observed = cache$observed + 1L else cache$fisher = cache$fisher + 1L
    dbreg:::dbglm_fe_step(conn, working, ..., raw_state = raw_state)
  }
  environment(driver) = e
  for (method in c("fisher", "newton")) {
    cache$observed = cache$fisher = 0L
    fit = driver(conn, "base", "x", "g", 4L, nrow(d), 0, list(rep(0, 4)),
      binomial("cloglog"), dbglm_control(irls = method))
    expect_true(fit$converged)
    expect_equal(fit$beta, unname(coef(oracle)[["x"]]), tolerance = 1e-8)
    expect_true(cache$fisher > 0L)
    if (method == "newton") expect_true(cache$observed > 0L) else expect_equal(cache$observed, 0L)
  }
  DBI::dbRemoveTable(conn, "base")
  expect_identical(DBI::dbListTables(conn), character())
})

# The identity includes the residual of the ORIGINAL group equations. It is
# exact even for deliberately inexact coefficients, independently of a solver.
local({
  set.seed(8832)
  F = cbind(1, c(0, 0, 1, 1, 2, 2), c(1, 0, 1, 0, 1, 0))
  Y = matrix(rnorm(18), 6, 3)
  w = c(.25, 1, 2, .5, 3, 1)
  H = crossprod(F, F * w)
  B = crossprod(F, Y * w)
  A = matrix(rnorm(9), 3, 3)
  information = dbreg:::fe_fwl_information
  e = new.env(parent = environment(information))
  e$fe_multiply = function(sys, coefficients) H %*% coefficients
  environment(information) = e
  global = list(cross = crossprod(Y, Y * w), n = 6, sum_w = sum(w))
  actual = information(NULL, "unused", c("v1", "v2", "v3"), B, A,
    list(sizes = 3L), NULL, dbreg:::projection_limits(), global = global)
  residual = Y - F %*% A
  expected = crossprod(residual, residual * w)
  expect_equal(unname(actual$within), unname(expected), tolerance = 1e-12)
  # Omitting the equation residual produces a materially wrong answer.
  wrong = global$cross - (crossprod(B, A) + crossprod(A, B)) / 2
  expect_true(max(abs(wrong - expected)) > 1)
})

# A well-scaled-looking result can still be small relative to its subtraction
# operands. The uncertainty guard must protect this case before rank selection.
local({
  calculate = dbreg:::dbglm_fe_information
  e = new.env(parent = environment(calculate))
  seen = new.env(parent = emptyenv())
  seen$recomputed = FALSE
  e$fe_information_query = function(...) {
    seen$recomputed = TRUE
    list(cross = diag(2), n = 4L)
  }
  e$sql_stable_sum = function(...) "SUM(w)"
  # The query uses a real tiny table; only the uncertain Gram is controlled.
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "within_rows", data.frame(w = rep(1, 4)))
  environment(calculate) = e
  H = matrix(c(1, 1 - 1e-7, 1 - 1e-7, 1), 2)
  raw = diag(c(1, 5e7, 5e7))
  M = rbind(c(1, 0, 0), cbind(0, H))
  result = calculate(conn, "within_rows", c("v1", "v2", "v3"),
    list(cross = raw, within = M),
    list2env(list(rows = 0L, cells = 0L), parent = emptyenv()),
    dbreg:::projection_limits(), 4L)
  expect_true(seen$recomputed)
  expect_equal(result$H, diag(2))
  DBI::dbRemoveTable(conn, "within_rows")
})

# Invalid normalization intermediates are uncertain information, so they must
# request actual-row recomputation before eigen(), including extreme units.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(w = rep(1, 4), v1 = c(-1, 1, -1, 1),
    v2 = c(-1, -1, 1, 1), v3 = c(-1, 1, 1, -1))
  DBI::dbWriteTable(conn, "within_rows", d)
  cases = list(
    list(H = matrix(c(5e-324, 1, 1, 5e-324), 2), raw = diag(c(1, 5e-324, 5e-324))),
    list(H = diag(rep(1e300, 2)), raw = diag(c(1, 1e308, 1e308))))
  for (case in cases) {
    M = rbind(c(1, 0, 0), cbind(0, case$H))
    value = dbreg:::dbglm_fe_information(conn, "within_rows", c("v1", "v2", "v3"),
      list(cross = case$raw, within = M),
      list2env(list(rows = 0L, cells = 0L), parent = emptyenv()),
      dbreg:::projection_limits(), 4L)
    expect_true(value$recomputed)
    expect_equal(unname(value$H), diag(rep(4, 2)), tolerance = 1e-12)
    expect_equal(value$score, c(0, 0), tolerance = 1e-12)
  }
  DBI::dbRemoveTable(conn, "within_rows")
  expect_identical(DBI::dbListTables(conn), character())
})
