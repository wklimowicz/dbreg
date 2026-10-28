# One diagonal FE family uses the existing packet and projection algebra.
# This is a correctness fixture, not a small-model performance benchmark.
local({
  library(DBI)
  info = list(list(id = 1:3, h = c(2, 5, 9)))
  system = dbreg:::fe_system(info, list())
  expect_identical(dim(system[["Cm"]]), c(3L, 0L))
  expect_identical(dim(system[["E"]]), c(0L, 0L))
  expect_identical(system[["rest"]], integer())
  solver = dbreg:::fe_solver(system, new.env(parent = emptyenv()), retain_null = TRUE)
  rhs = matrix(c(1, 3, -2, 4, 7, 2), 3L, 2L)
  expect_equal(solver[["solve"]](rhs), rhs / c(2, 5, 9), tolerance = 1e-14)
  expect_identical(dim(solver[["solve"]](matrix(0, 3L, 0L))), c(3L, 0L))
  expect_identical(solver[["project"]](rhs), rhs)
  expect_null(solver[["null_basis"]])

  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")), bigint = "integer64")
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  p = data.frame(w = c(.5, 2, 1, 3, .75, 1.25, 4, 2), x = c(1, 3, -2, 4, 0, 7, 2, -1),
    q = c(2, -1, 3, 0, 4, 2, -3, 5), z = c(1, -2, 4, 3, 0, 1, 2, -1), score = c(.2, -.1, 1, -2, .5, 0, 3, -.4))
  group = c(2L, 1L, 2L, 3L, 1L, 2L, 3L, 3L)
  keys = bit64::as.integer64(c("9007199254740992", "9007199254740993", "9007199254740994"))
  p[["g"]] = keys[group]
  DBI::dbWriteTable(conn, "k1_packet", p)
  cols = c("x", "q", "z")
  packet = dbreg:::query_fe_moments(conn, "k1_packet", "g", cols, extra_sums = "score", global_terms = "COUNT(*) AS n_rows")
  expect_identical(as.character(packet[["info"]][[1L]][["id"]]), as.character(keys))
  expect_equal(packet[["rows"]], 3)
  expect_equal(as.numeric(packet[["globals"]][["n_rows"]]), nrow(p))
  expect_identical(dim(packet[["sys"]][["Cm"]]), c(3L, 0L))
  expect_identical(dim(packet[["sys"]][["E"]]), c(0L, 0L))
  raw = as.matrix(p[cols])
  group_w = drop(rowsum(p[["w"]], group, reorder = TRUE))
  means = rowsum(p[["w"]] * raw, group, reorder = TRUE) / group_w
  centered = raw - means[group, , drop = FALSE]
  expect_equal(unname(packet[["B"]]), unname(rowsum(p[["w"]] * raw, group, reorder = TRUE)), tolerance = 1e-13)
  expect_equal(as.numeric(packet[["B_extra"]]), as.numeric(rowsum(p[["score"]], group, reorder = TRUE)), tolerance = 1e-13)
  control = dbreg::dbglm_control()
  forbid = function(...) stop("Unexpected physical projection request in a moment-only test.")
  projected = dbreg:::project_fe_moments(conn, "k1_packet", cols, "g", list(list(f = "g", v = "1.0")),
    sqrt(colSums(p[["w"]] * raw^2)), control[["fe_tol"]], control, forbid, forbid,
    new.env(parent = emptyenv()), materialize = FALSE, extra_sums = "score", need_table = FALSE, precomputed = packet)
  expect_equal(unname(projected[["coefficients"]]), unname(means), tolerance = 1e-13)
  expect_equal(unname(crossprod(raw, p[["w"]] * raw) - projected[["cross_reduction"]]),
    unname(crossprod(centered, p[["w"]] * centered)), tolerance = 1e-12)
  expect_equal(as.numeric(crossprod(raw, p[["score"]]) - projected[["extra_reduction"]]),
    as.numeric(crossprod(centered, p[["score"]])), tolerance = 1e-12)
  empty_packet = dbreg:::query_fe_moments(conn, "k1_packet", "g", character())
  expect_identical(dim(empty_packet[["B"]]), c(3L, 0L))
  expect_identical(dim(empty_packet[["B_extra"]]), c(3L, 0L))
  expect_identical(dim(dbreg:::fe_solver(empty_packet[["sys"]], new.env(parent = emptyenv()))[["solve"]](empty_packet[["B"]])), c(3L, 0L))
  wire = dbreg:::read_fe_moments(conn, dbreg:::fe_moments_sql(conn, "k1_packet", "g", dbreg:::fe_value_terms(conn, cols, "score"), "COUNT(*) AS n_rows", compact = TRUE))
  decoded = dbreg:::decode_fe_moments(wire, 1L, 3L, 1L)
  expect_identical(decoded[["pairs"]], list())
  expect_equal(sort(as.integer(wire[["__set"]])), c(0L, 1L, 1L, 1L))

  # Native K1 must decline before query/cache enablement, including a stale
  # enabled cache. Reuse the existing local native fixture; its subclass makes
  # any SQL boundary a hard failure rather than silently accepting Duck syntax.
  source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
  remote = dbglm_test_databricks_connection(conn)
  no_query_class = "dbreg_k1_no_query_fixture"
  if (!methods::isClass(no_query_class)) methods::setClass(no_query_class, contains = "dbglm_sql_fixture")
  calls = new.env(parent = emptyenv())
  calls[["sql"]] = character()
  no_sql = function(conn, statement, ...) {
    calls[["sql"]] = c(calls[["sql"]], statement)
    stop("Native single-FE eligibility must decline before SQL.")
  }
  methods::setMethod("dbGetQuery", c(no_query_class, "character"), no_sql)
  methods::setMethod("dbExecute", c(no_query_class, "character"), no_sql)
  no_query = methods::new(no_query_class, native = conn)
  for (enabled in c(FALSE, TRUE)) {
    cache = new.env(parent = emptyenv())
    if (enabled) cache[["moment_enabled"]] = TRUE
    expect_false(dbreg:::dbglm_moments_enabled(no_query, "missing_source", "g", control, cache, 3L))
    expect_identical(cache[["moment_enabled"]], FALSE)
    expect_identical(calls[["sql"]], character())
  }

  set.seed(281028)
  n = 480L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:8, each = 60L),
    cl = sample(1:16, n, TRUE), w = sample(1:3, n, TRUE), o = runif(n, -.1, .1), t = rnorm(n))
  eta = -.3 + .4 * d[["x"]] - .2 * d[["z"]] + sin(d[["g"]]) / 5 + d[["o"]]
  d[["y"]] = rbinom(n, 1L, plogis(eta))
  d[["yp"]] = rpois(n, exp(eta / 2))
  d[["yg"]] = eta + rnorm(n) / 2
  DBI::dbWriteTable(conn, "k1_values", d)
  inventory = sort(DBI::dbListTables(conn))
  reference = stats::glm(y ~ x + z + factor(g), data = d, weights = w, offset = o,
    family = stats::binomial(), control = stats::glm.control(epsilon = 1e-12, maxit = 100L))
  X = stats::model.matrix(reference)
  mu = stats::fitted(reference)
  bread = solve(crossprod(X, (d[["w"]] * mu * (1 - mu)) * X))
  scores = rowsum((d[["w"]] * (d[["y"]] - mu)) * X, d[["cl"]], reorder = FALSE)
  G = nrow(scores); K = ncol(X)
  V = bread %*% crossprod(scores) %*% bread * G / (G - 1) * (n - 1) / (n - K)
  close_check = function(fit) {
    close(fit)
    expect_true(DBI::dbIsValid(conn))
    expect_identical(sort(DBI::dbListTables(conn)), inventory)
  }
  fit_case = function(control = dbreg::dbglm_control(), connection = conn) dbreg::dbglm(y ~ x + z | g,
    conn = connection, table = "k1_values", family = stats::binomial(), weights = "w", offset = "o", vcov = ~cl,
    fixef_rm = "none", strategy = "moments", retain = "prediction", control = control)
  check_logit = function(fit) {
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(reference)[c("x", "z")]), tolerance = 1e-7)
    expect_identical(dim(vcov(fit)), c(2L, 2L))
    expect_equal(as.numeric(vcov(fit)), as.numeric(V[c("x", "z"), c("x", "z")]), tolerance = 1e-7)
    expect_equal(nobs(fit), n)
    expect_equal(fit[["n_params_model"]], K)
    expect_equal(fit[["n_params_vcov"]], K)
    expect_equal(stats::df.residual(fit), n - K)
    expect_equal(as.numeric(logLik(fit)), sum(d[["w"]] * stats::dbinom(d[["y"]], 1L, mu, log = TRUE)), tolerance = 1e-8)
    expect_equal(fit[["deviance"]], stats::deviance(reference), tolerance = 1e-8)
    expect_equal(as.numeric(predict(fit, newdata = d, type = "response")), as.numeric(mu), tolerance = 1e-7)
    expect_true(fit[["score_max"]] <= fit[["control"]][["score_tol"]])
    expect_true(fit[["ap_residual"]] <= fit[["control"]][["fe_tol"]])
    expect_true(all(diff(fit[["deviance_history"]]) <= 1e-9))
    expect_true(fit[["step_halvings"]] >= 0L)
  }
  fast = fit_case()
  expect_identical(fast[["projection_method"]], "moments")
  check_logit(fast); close_check(fast)
  native_group = fit_case(connection = remote)
  expect_identical(native_group[["projection_method"]], "group")
  check_logit(native_group); close_check(native_group)
  sql = fit_case(dbreg::dbglm_control(projection_method = "sql"))
  expect_identical(sql[["projection_method"]], "group")
  check_logit(sql); close_check(sql)
  declined = local({
    saved = options(dbreg.duckdb_ap_max_entries = 1)
    on.exit(options(saved))
    fit_case()
  })
  expect_identical(declined[["projection_method"]], "group")
  expect_true(any(vapply(declined[["projection_diagnostics"]], function(attempt)
    identical(attempt[["status"]], "declined") && identical(attempt[["reason"]], "system_limit"), logical(1))))
  check_logit(declined); close_check(declined)

  # Existing independent family reference; no package moments enter the oracle.
  source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
  for (family_case in list(list(y = "yp", family = stats::poisson(), vcov = "hc1", ref_vcov = "hetero"),
    list(y = "yg", family = stats::gaussian(), vcov = "iid", ref_vcov = "iid"))) {
    formula = stats::as.formula(paste(family_case[["y"]], "~x+z|g"))
    fit = dbreg::dbglm(formula, conn = conn, table = "k1_values", family = family_case[["family"]],
      weights = "w", offset = "o", vcov = family_case[["vcov"]], fixef_rm = "none", strategy = "moments")
    ref = dbglm_reference_feglm(formula, d, family = family_case[["family"]], weights = ~w, offset = ~o,
      vcov = family_case[["ref_vcov"]], fixef.rm = "none")
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-7)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-8)
    close_check(fit)
  }
  only = dbreg::dbglm(y ~ 1 | g, conn = conn, table = "k1_values", family = stats::binomial(),
    weights = "w", vcov = "iid", fixef_rm = "none", retain = "prediction", strategy = "moments")
  group_mu = drop(rowsum(d[["w"]] * d[["y"]], d[["g"]])) / drop(rowsum(d[["w"]], d[["g"]]))
  expect_equal(length(coef(only)), 0L)
  expect_identical(dim(vcov(only)), c(0L, 0L))
  expect_equal(as.numeric(predict(only, newdata = d, type = "response")), as.numeric(group_mu[d[["g"]]]), tolerance = 1e-7)
  close_check(only)
  no_intercept = dbreg::dbglm(y ~ 0 + x + z | g, conn = conn, table = "k1_values", family = stats::binomial(),
    weights = "w", offset = "o", vcov = ~cl, fixef_rm = "none", strategy = "moments")
  expect_equal(as.numeric(coef(no_intercept)), as.numeric(coef(reference)[c("x", "z")]), tolerance = 1e-7)
  expect_equal(as.numeric(vcov(no_intercept)), as.numeric(V[c("x", "z"), c("x", "z")]), tolerance = 1e-7)
  close_check(no_intercept)
  trend = dbreg::dbglm(yg ~ x | g[t], conn = conn, table = "k1_values", family = stats::gaussian(),
    weights = "w", offset = "o", fixef_rm = "none", strategy = "moments")
  trend_ref = dbglm_reference_feglm(yg ~ x | g[t], d, family = stats::gaussian(), weights = ~w, offset = ~o, fixef.rm = "none")
  expect_identical(trend[["projection_method"]], "group")
  expect_equal(as.numeric(coef(trend)), as.numeric(coef(trend_ref)), tolerance = 1e-7)
  expect_equal(as.numeric(vcov(trend)), as.numeric(vcov(trend_ref)), tolerance = 1e-7)
  close_check(trend)

  cells = expand.grid(x = c(-1, 0, 1), z = c(-.5, .5), g = 1:4, o = c(-.1, .1))
  cd = cells[rep(seq_len(nrow(cells)), each = 12L), , drop = FALSE]
  cd[["w"]] = rep(1:3, length.out = nrow(cd))
  cd[["cl"]] = rep(1:4, length.out = nrow(cd))
  cd[["y"]] = rbinom(nrow(cd), 1L, plogis(-.3 + .4 * cd[["x"]] - .2 * cd[["z"]] + sin(cd[["g"]]) / 5 + cd[["o"]]))
  compressed = dbreg::dbglm(y ~ x + z | g, conn = conn, data = cd, family = stats::binomial(),
    weights = "w", offset = "o", vcov = ~cl, fixef_rm = "none", strategy = "compress", retain = "prediction")
  compressed_ref = dbglm_reference_feglm(y ~ x + z | g, cd, family = stats::binomial(), weights = ~w, offset = ~o,
    vcov = ~cl, fixef.rm = "none")
  expect_true(compressed[["n_cells"]] < nrow(cd))
  expect_equal(as.numeric(coef(compressed)), as.numeric(coef(compressed_ref)), tolerance = 1e-7)
  expect_equal(as.numeric(vcov(compressed)), as.numeric(vcov(compressed_ref)), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(compressed)), as.numeric(logLik(compressed_ref)), tolerance = 1e-8)
  expect_equal(as.numeric(predict(compressed, newdata = cd, type = "response")), as.numeric(predict(compressed_ref, newdata = cd, type = "response")), tolerance = 1e-7)
  close_check(compressed)
})
