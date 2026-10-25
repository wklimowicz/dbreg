library(DBI)
library(dbreg)
library(fixest)

(function() {
  old = options(dbreg.duckdb_ap_method = "cg", dbreg.ap_tol = 1e-10,
    dbreg.ap_max_iter = 100L, dbreg.duckdb_ap_max_edges = 1e6,
    dbreg.duckdb_ap_max_entries = 2e6)
  con = dbConnect(duckdb::duckdb(),
    config = list(temp_directory = "", memory_limit = "1GiB", threads = "2"))
  on.exit({ options(old); dbDisconnect(con, shutdown = TRUE) })
  set.seed(20261022)
  n = 1800L
  dat = data.frame(fe1 = sample(1:53, n, TRUE), fe2 = sample(1:29, n, TRUE),
    fe3 = sample(1:11, n, TRUE), cluster = sample(1:37, n, TRUE),
    x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 0.1, 3),
    category = sample(c("control", "other", "treated"), n, TRUE))
  dat$y = .7 * dat$x1 - .4 * dat$x2 + .2 * dat$x1 * dat$x2 +
    sin(dat$fe1) + cos(dat$fe2) + .1 * sin(dat$fe3) + rnorm(n, sd = .3)
  dbWriteTable(con, "cg_fixture", dat)
  no_leaks = function() expect_false(any(startsWith(dbListTables(con), "dbreg_ap_")))
  check = function(fml, weighted = FALSE, covariance = "iid") {
    args = list(fml = fml, conn = con, table = "cg_fixture", strategy = "demean",
      vcov = covariance, ssc = "full")
    ref_args = list(fml = fml, data = dat, fixef.rm = "none",
      ssc = fixest::ssc(K.fixef = "full"),
      vcov = if (identical(covariance, "hc1")) "hetero" else covariance)
    if (weighted) {
      args$weights = "w"
      ref_args$weights = ~w
    }
    got = do.call(dbreg, args)
    ref = do.call(fixest::feols, ref_args)
    ref_names = sub("^I\\(I\\((.*)\\)\\)$", "I(\\1)", names(coef(ref)))
    order = match(names(coef(got)), ref_names)
    expect_true(!anyNA(order))
    expect_true(max(abs(as.numeric(coef(got)) - as.numeric(coef(ref))[order])) < 1e-6)
    expect_true(max(abs(vcov(got) - unclass(vcov(ref))[order, order])) < 1e-6)
    expect_identical(got$ap_method, "bounded_pcg")
    expect_true(got$ap_max_abs <= 1e-10)
    expect_true(got$ap_iterations <= 100L)
    no_leaks()
    got
  }
  simple = y ~ x1 + x2 | fe1 + fe2
  for (weighted in c(FALSE, TRUE)) {
    for (covariance in list("iid", "hc1", ~cluster)) check(simple, weighted, covariance)
  }
  check(y ~ x1*x2 + category | fe1 + fe2, TRUE, ~cluster)

  # Distinct BIGINT FE keys may round to the same R double. Projection keys
  # must retain all 53 groups and produce the same estimates as integer FE.
  dbExecute(con, paste(
    "CREATE TABLE cg_big AS SELECT *,",
    "9007199254740992::BIGINT + fe1::BIGINT AS bigfe FROM cg_fixture"
  ))
  big = dbreg(y ~ x1 + x2 | bigfe + fe2, conn = con, table = "cg_big",
    strategy = "demean", weights = "w", vcov = ~cluster)
  standard = check(simple, TRUE, ~cluster)
  expect_identical(big$ap_method, "bounded_pcg")
  expect_equal(coef(big), coef(standard), tolerance = 1e-8)
  expect_true(max(abs(vcov(big) - vcov(standard))) < 1e-8)
  expect_equal(big$nobs_orig, n)
  expect_equal(big$n_fe_levels, c(53L, 29L))
  no_leaks()

  # String FE and disconnected graphs remain valid. Compare with the join
  # implementation's covariance convention for disconnected FE dimensions.
  dbExecute(con, paste(
    "CREATE TABLE cg_disconnected AS SELECT * EXCLUDE(fe2),",
    "CAST(fe2 + CASE WHEN fe1 <= 26 THEN 0 ELSE 100 END AS VARCHAR) AS fe2",
    "FROM cg_fixture"
  ))
  fit_method = function(table, method) {
    options(dbreg.duckdb_ap_method = method)
    dbreg(simple, conn = con, table = table, weights = "w", vcov = ~cluster,
      strategy = "demean")
  }
  fast = fit_method("cg_disconnected", "cg")
  joined = fit_method("cg_disconnected", "join")
  expect_identical(fast$ap_method, "bounded_pcg")
  expect_equal(coef(fast), coef(joined), tolerance = 1e-7)
  expect_true(max(abs(vcov(fast) - vcov(joined))) < 1e-7)
  no_leaks()

  # Transfer limits are checked before returning FE summaries to R, and
  # rejected attempts clean up their state before the existing join fallback.
  for (caps in list(list(dbreg.duckdb_ap_max_edges = 1L),
                    list(dbreg.duckdb_ap_max_entries = 1L))) {
    options(dbreg.duckdb_ap_max_edges = 1e6, dbreg.duckdb_ap_max_entries = 2e6)
    options(caps)
    bounded = fit_method("cg_fixture", "cg")
    expect_identical(bounded$ap_method, "row_tables")
    expect_equal(coef(bounded), coef(standard), tolerance = 1e-8)
    no_leaks()
  }
  options(dbreg.duckdb_ap_max_edges = 1e6, dbreg.duckdb_ap_max_entries = 2e6)
  expect_error({
    options(dbreg.duckdb_ap_max_edges = 1.5)
    fit_method("cg_fixture", "cg")
  }, "positive whole numbers")
  options(dbreg.duckdb_ap_max_edges = 1e6)

  # Nonconvergence and SQL failures must not leave temporary state behind.
  expect_null(dbreg:::dbreg_duckdb_cg_projections(con, "FROM cg_fixture",
    c("fe1", "fe2"), "y", c("x1", "x2"), c("x1", "x2"), "w",
    max_iter = 1L, tol = 1e-14))
  no_leaks()
  expect_error(dbreg:::dbreg_duckdb_cg_projections(con, "FROM cg_fixture",
    c("fe1", "fe2"), "y", "missing_column", "x1", NULL), "missing_column")
  no_leaks()

  # More than two FE use the established implementation with the option on.
  options(dbreg.duckdb_ap_method = "cg")
  multi = dbreg(y ~ x1 + x2 | fe1 + fe2 + fe3, conn = con,
    table = "cg_fixture", weights = "w", strategy = "demean")
  ref = fixest::feols(y ~ x1 + x2 | fe1 + fe2 + fe3, data = dat, weights = ~w,
    fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
  expect_identical(multi$ap_method, "row_tables")
  expect_equal(coef(multi), coef(ref), tolerance = 1e-6)
  expect_true(max(abs(vcov(multi) - unclass(vcov(ref)))) < 1e-6)
  no_leaks()
  options(dbreg.duckdb_ap_method = "invalid")
  expect_error(dbreg(simple, conn = con, table = "cg_fixture", strategy = "demean"),
    'must be "join" or "cg"')
  no_leaks()
})()
