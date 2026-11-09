library(DBI)

# A projection must preserve the raw score meat, including IV's distinct
# score coordinates, multiple residuals and the intercept-only case.
local({
  if (!methods::isClass("dbreg_hc1_sql_fixture")) {
    methods::setClass("dbreg_hc1_sql_fixture", contains = "DBIConnection", slots = c(native = "ANY", queries = "environment"))
    methods::setMethod("dbGetInfo", "dbreg_hc1_sql_fixture", function(dbObj, ...) list(dbms.name = "duckdb"))
    methods::setMethod("dbGetQuery", c("dbreg_hc1_sql_fixture", "character"), function(conn, statement, ...) {
      queries = methods::slot(conn, "queries")
      queries[["sql"]] = c(queries[["sql"]], statement)
      DBI::dbGetQuery(methods::slot(conn, "native"), statement, ...)
    })
  }
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  queries = new.env(parent = emptyenv())
  conn = methods::new("dbreg_hc1_sql_fixture", native = native, queries = queries)
  set.seed(1809)
  d = data.frame(y = rnorm(120), x = rnorm(120), z = rnorm(120), w = runif(120, .2, 2), g = rep(1:12, each = 10))
  # Source names deliberately overlap the generated projection aliases.
  d[["__dbreg_hc1_r1"]] = 999
  d[["__dbreg_hc1_s2"]] = 999
  d[["__dbreg_hc1_w_sq"]] = 999
  DBI::dbWriteTable(native, "hc1_source", d)
  inputs = list2env(list(conn = conn, cluster_var = "g", vcov_type_req = "hc1"))
  frame = list(cte_sql = "WITH base AS (SELECT * FROM hc1_source)", cte_name = "base", weights_expr = NULL)
  residuals = c(model = "(y - (0.7 + 0.2 * x))", first_stage = "(x - (0.1 + 0.4 * z))")
  errors = list(d[["y"]] - (.7 + .2*d[["x"]]), d[["x"]] - (.1 + .4*d[["z"]]))
  for (weighted in c(FALSE, TRUE)) {
    frame[["weights_expr"]] = if (weighted) "w" else NULL
    w = if (weighted) d[["w"]] else rep(1, nrow(d))
    for (scores in list(character(), "z", c("x", "z"))) {
      X = cbind(1, as.matrix(d[, scores, drop = FALSE]))
      colnames(X) = c("(Intercept)", scores)
      before = length(queries[["sql"]])
      result = dbreg:::compute_score_meat(inputs, frame, residuals, scores, score_intercept = TRUE)
      expect_identical(length(queries[["sql"]]) - before, 1L)
      expect_identical(names(result), names(residuals))
      for (r in seq_along(errors)) expect_equal(result[[r]], crossprod(X * (errors[[r]] * w)), tolerance = 1e-11)
    }
    # Unchanged cluster path: an independent grouped-score calculation.
    inputs[["vcov_type_req"]] = "cluster"
    result = dbreg:::compute_score_meat(inputs, frame, residuals, "z", score_intercept = TRUE)
    X = cbind(1, z = d[["z"]]); colnames(X)[1] = "(Intercept)"
    for (r in seq_along(errors)) {
      expected = crossprod(rowsum(X * (errors[[r]] * w), d[["g"]]))
      attr(expected, "n_clusters") = 12L
      expect_equal(result[[r]], expected, tolerance = 1e-11)
    }
    inputs[["vcov_type_req"]] = "hc1"
  }
  # Squaring residuals before multiplying scores is intentional: forming
  # row scores first would change underflow/overflow at extreme scales.
  for (exponent in c(-200, 200)) {
    frame = list(cte_sql = sprintf("WITH base AS (SELECT 1e%d AS e, 1e%d AS x)", exponent, -exponent), cte_name = "base", weights_expr = NULL)
    result = dbreg:::compute_score_meat(inputs, frame, c(model = "e"), "x")
    expect_identical(as.numeric(result[[1]]), if (exponent < 0) 0 else Inf)
  }
  # A wide residual is written once, rather than repeated in every matrix
  # entry. This is a SQL-size bound, not an elapsed-time assertion.
  wide = as.data.frame(matrix(rnorm(1600), 50, 32))
  names(wide) = paste0("x", seq_len(32))
  wide[["y"]] = rnorm(50)
  DBI::dbWriteTable(native, "hc1_wide", wide)
  frame = list(cte_sql = "WITH base AS (SELECT * FROM hc1_wide)", cte_name = "base", weights_expr = NULL)
  residual = paste0("(y - (", paste(sprintf("0.0123456789012345 * x%d", 1:32), collapse = " + "), "))")
  result = dbreg:::compute_score_meat(inputs, frame, c(model = residual), names(wide)[1:32], score_intercept = TRUE)
  expect_true(nchar(tail(queries[["sql"]], 1), type = "bytes") < 100000L)
  X = cbind(1, as.matrix(wide[, 1:32]))
  e = wide[["y"]] - rowSums(as.matrix(wide[, 1:32])) * .0123456789012345
  expect_equal(as.numeric(result[[1]]), as.numeric(crossprod(X * e)), tolerance = 1e-10)
  expect_identical(sort(DBI::dbListTables(native)), c("hc1_source", "hc1_wide"))
})

# Public weighted HC1 regression/IV contracts against independent estimators.
local({
  set.seed(921)
  d = data.frame(x = rnorm(180), z = rnorm(180), w = runif(180, .4, 2))
  d[["endog"]] = .4*d[["x"]] + d[["z"]] + rnorm(180)
  d[["y"]] = .7 + .4*d[["x"]] + .9*d[["endog"]] + rnorm(180)
  for (formula in list(y ~ x, y ~ x | endog ~ z)) {
    fit = if (length(formula) == 3 && !grepl("|", deparse(formula), fixed = TRUE)) dbreg else dbivreg
    actual = fit(formula, data = d, weights = "w", vcov = "hc1", strategy = "moments")
    expected = fixest::feols(formula, data = d, weights = ~w, vcov = "hc1")
    b = coef(expected)
    names(b) = sub("^fit_", "", names(b))
    V = vcov(expected)
    dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
    terms = names(coef(actual))
    expect_equal(as.numeric(coef(actual)), as.numeric(b[terms]), tolerance = 1e-9)
    expect_equal(as.numeric(vcov(actual)), as.numeric(V[terms, terms]), tolerance = 1e-9)
  }
})
