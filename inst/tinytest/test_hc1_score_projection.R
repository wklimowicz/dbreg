# Independent score/sandwich references and an emitted-work bound, not timing.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  set.seed(14192)
  n = 180L
  d = data.frame(x = rnorm(n), z = rnorm(n), u = rnorm(n),
    w = runif(n, .3, 1.7), o = runif(n, -.2, .2),
    g = rep(seq_len(12L), 15L), h = rep(seq_len(9L), each = 20L))
  d[["instrument"]] = rnorm(n)
  d[["endog"]] = .8 * d[["instrument"]] + .3 * d[["x"]] + rnorm(n)
  d[["y"]] = .4 + .7 * d[["x"]] - .2 * d[["z"]] + .5 * d[["endog"]] +
    .1 * sin(d[["g"]]) + d[["o"]] + d[["u"]]
  # These physical names must remain usable inside the private projection.
  d[["__dbreg_hc1_r1"]] = d[["z"]]
  d[["__dbreg_hc1_s1"]] = d[["x"]]
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "hc1_score_source", d)
  inputs = function(conn, type = "hc1") list2env(list(conn = conn,
    vcov_type_req = type, cluster_var = "g"), parent = emptyenv())
  frame = function(weighted, scale = 1) list(
    cte_sql = "WITH score_base AS (SELECT * FROM hc1_score_source)",
    cte_name = "score_base", weights_expr = if (weighted) "w" else NULL,
    weight_scale = scale)
  residuals = c(first = "y-0.4-0.7*x+0.2*z", second = "u+0.3*z")
  E = cbind(d[["y"]] - .4 - .7 * d[["x"]] + .2 * d[["z"]], d[["u"]] + .3 * d[["z"]])
  S = cbind(1, d[["x"]], d[["z"]])
  score_names = c("x", "z")
  for (backend in list(conn, dbglm_test_databricks_connection(conn))) {
    for (weighted in c(FALSE, TRUE)) {
      scale = if (weighted) max(d[["w"]]) else 1
      a = if (weighted) d[["w"]] / scale else rep(1, n)
      value = dbreg:::compute_score_meat(inputs(backend), frame(weighted, scale),
        residuals, c("x", "z"), score_intercept = TRUE)
      expect_identical(names(value), names(residuals))
      for (r in seq_len(ncol(E))) {
        reference = crossprod(S, S * (a * E[, r])^2)
        expect_equal(as.numeric(value[[r]]), as.numeric(reference), tolerance = 1e-12)
        expect_identical(dimnames(value[[r]]), rep(list(c("(Intercept)", score_names)), 2L))
        expect_equal(attr(value[[r]], "weight_scale"), scale)
      }
      # A constant-only score omits both score factors, just as an intercept did.
      only = dbreg:::compute_score_meat(inputs(backend), frame(weighted, scale),
        residuals, character(), score_intercept = TRUE)
      expect_equal(as.numeric(only[[1L]]), sum((a * E[, 1L])^2), tolerance = 1e-12)
      expect_identical(dim(only[[1L]]), c(1L, 1L))
      # The cluster path is an unchanged control, independently grouped in R.
      cluster = dbreg:::compute_score_meat(inputs(backend, "cluster"), frame(weighted, scale),
        residuals, c("x", "z"), score_intercept = TRUE)
      for (r in seq_len(ncol(E))) {
        reference = crossprod(rowsum(S * (a * E[, r]), d[["g"]]))
        expect_equal(as.numeric(cluster[[r]]), as.numeric(reference), tolerance = 1e-12)
        expect_equal(attr(cluster[[r]], "n_clusters"), length(unique(d[["g"]])))
      }
    }
    collision = dbreg:::compute_score_meat(inputs(backend), frame(FALSE),
      c(first = "y-0.4-0.7*x+0.2*z"), c("__dbreg_hc1_s1", "__dbreg_hc1_r1"))
    expected = crossprod(S[, -1L], S[, -1L] * E[, 1L]^2)
    expect_equal(as.numeric(collision[[1L]]), as.numeric(expected), tolerance = 1e-12)
  }
  # Low-level extreme factors distinguish the specified left-to-right order
  # from pre-multiplied scores, which would turn both results into one.
  scalar = list(cte_sql = "WITH scalar_source AS (SELECT 1 AS marker)",
    cte_name = "scalar_source", weights_expr = NULL, weight_scale = 1)
  over = dbreg:::compute_score_meat(inputs(conn), scalar, c(r = "1e200"), "1e-200")
  under = dbreg:::compute_score_meat(inputs(conn), scalar, c(r = "1e-200"), "1e200")
  expect_true(is.infinite(over[[1L]][1L, 1L]) && over[[1L]][1L, 1L] > 0)
  expect_identical(under[[1L]][1L, 1L], 0)
  # Normalize the original weight BEFORE squaring, including at extreme scale.
  for (magnitude in c(1e-150, 1e150, 1e300)) {
    scale = max(d[["w"]]) * magnitude
    ext = frame(TRUE, scale)
    ext[["weights_expr"]] = paste0("(w*", dbreg:::dbglm_number(magnitude), ")")
    actual = dbreg:::compute_score_meat(inputs(conn), ext, residuals, c("x", "z"), score_intercept = TRUE)
    reference = crossprod(S, S * ((d[["w"]] / max(d[["w"]])) * E[, 1L])^2)
    expect_equal(as.numeric(actual[[1L]]), as.numeric(reference), tolerance = 1e-12)
  }
  # Public OLS: full QR and HC1 sandwich, with offset and explicit SSC control.
  for (weighted in c(FALSE, TRUE)) for (adjust in c(FALSE, TRUE)) {
    a = if (weighted) d[["w"]] else rep(1, n)
    beta = qr.coef(qr(S * sqrt(a)), (d[["y"]] - d[["o"]]) * sqrt(a))
    e = as.numeric(d[["y"]] - d[["o"]] - S %*% beta)
    bread = solve(crossprod(S, S * a))
    V = bread %*% crossprod(S, S * (a * e)^2) %*% bread
    if (adjust) V = V * n / (n - ncol(S))
    fit = dbreg(y ~ x + z, data = d, weights = if (weighted) "w" else NULL,
      offset = "o", vcov = "hc1", strategy = "moments", ssc = list(K.adj = adjust))
    expect_equal(as.numeric(coef(fit)), as.numeric(beta), tolerance = 1e-11)
    expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-11)
    expect_equal(fit[["n_params"]], ncol(S))
    expect_equal(fit[["df_t"]], n - ncol(S))
    close(fit)
  }
  # Absorbed weighted FE: a full dummy QR is independent of the SQL projection.
  X = model.matrix(~x + z + factor(g) + factor(h), d)
  a = d[["w"]]
  beta = qr.coef(qr(X * sqrt(a)), (d[["y"]] - d[["o"]]) * sqrt(a))
  e = as.numeric(d[["y"]] - d[["o"]] - X %*% beta)
  bread = solve(crossprod(X, X * a))
  V = n / (n - ncol(X)) * bread %*% crossprod(X, X * (a * e)^2) %*% bread
  fit = dbreg(y ~ x + z | g + h, data = d, weights = "w", offset = "o",
    vcov = "hc1", strategy = "demean")
  expect_equal(as.numeric(coef(fit)), as.numeric(beta[c("x", "z")]), tolerance = 1e-9)
  expect_equal(as.numeric(vcov(fit)), as.numeric(V[c("x", "z"), c("x", "z")]), tolerance = 1e-9)
  expect_equal(fit[["n_params"]], ncol(X))
  close(fit)
  # IV structural covariance: residual regressors differ from instrument scores.
  X = cbind(1, d[["x"]], d[["endog"]])
  Z = cbind(1, d[["x"]], d[["instrument"]], d[["z"]])
  ZZ = crossprod(Z, Z * a)
  XZ = crossprod(X, Z * a)
  projection = XZ %*% solve(ZZ)
  bread = solve(projection %*% t(XZ))
  beta = bread %*% projection %*% crossprod(Z, a * (d[["y"]] - d[["o"]]))
  e = as.numeric(d[["y"]] - d[["o"]] - X %*% beta)
  V = n / (n - ncol(X)) * bread %*% projection %*% crossprod(Z, Z * (a * e)^2) %*% t(projection) %*% bread
  fit = dbivreg(y ~ x | endog ~ instrument + z, data = d, weights = "w", offset = "o",
    vcov = "hc1", strategy = "moments")
  expect_equal(as.numeric(coef(fit)), as.numeric(beta), tolerance = 1e-9)
  expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-9)
  expect_equal(fit[["n_params"]], ncol(X))
  close(fit)
})

local({
  # Observe the actual S4 query method and execute it, without namespace tracing.
  if (!methods::isClass("dbreg_hc1_collector")) methods::setClass("dbreg_hc1_collector",
    contains = "DBIConnection", slots = c(native = "ANY", events = "environment"))
  methods::setMethod("dbGetInfo", "dbreg_hc1_collector", function(dbObj, ...) list(dbms.name = "DuckDB"))
  methods::setMethod("dbIsValid", "dbreg_hc1_collector", function(dbObj, ...) DBI::dbIsValid(methods::slot(dbObj, "native")))
  methods::setMethod("dbGetQuery", c("dbreg_hc1_collector", "character"), function(conn, statement, ...) {
    events = methods::slot(conn, "events")
    events[["queries"]] = c(events[["queries"]], list(statement))
    DBI::dbGetQuery(methods::slot(conn, "native"), statement, ...)
  })
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", temp_directory = "")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE), add = TRUE)
  events = new.env(parent = emptyenv())
  events[["queries"]] = list()
  conn = methods::new("dbreg_hc1_collector", native = native, events = events)
  set.seed(2904)
  X = matrix(rnorm(64L * 32L), 64L, 32L, dimnames = list(NULL, paste0("x", seq_len(32L))))
  d = as.data.frame(X)
  d[["y"]] = rowSums(X) + rnorm(nrow(X))
  DBI::dbWriteTable(native, "hc1_wide_source", d)
  before = sort(DBI::dbListTables(native))
  frame = list(cte_sql = "WITH wide_base AS (SELECT * FROM hc1_wide_source)", cte_name = "wide_base",
    weights_expr = NULL, weight_scale = 1)
  inputs = list2env(list(conn = conn, cluster_var = NULL, vcov_type_req = "hc1"), parent = emptyenv())
  residual = paste0("y-(", paste(colnames(X), collapse = "+"), ")")
  value = dbreg:::compute_score_meat(inputs, frame, c(model = residual), colnames(X), score_intercept = TRUE)
  S = cbind(1, X)
  reference = crossprod(S, S * (d[["y"]] - rowSums(X))^2)
  expect_equal(as.numeric(value[[1L]]), as.numeric(reference), tolerance = 1e-12)
  expect_equal(length(events[["queries"]]), 1L)
  sql = events[["queries"]][[1L]]
  # A finite emitted-work bound exposes accidental repetition of the predictor;
  # no elapsed-time threshold, exact SQL spelling or particular CTE is required.
  expect_true(nchar(sql, type = "bytes") < 100000L)
  expect_false(grepl("\\b(MATERIALIZED|CREATE|INSERT|DROP)\\b", sql, ignore.case = TRUE, perl = TRUE))
  expect_identical(sort(DBI::dbListTables(native)), before)
})
