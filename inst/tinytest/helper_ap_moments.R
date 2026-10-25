# Shared numerical/ownership checks for the bounded OLS/IV within transform.
dbreg_test_ap_moments = function(conn, verbose = FALSE) {
  saved = options(dbreg.ap_solver = "moments", dbreg.ap_tol = 1e-10,
    dbreg.ap_max_iter = 100L, dbreg.ap_method = NULL, dbreg.duckdb_ap_method = NULL)
  on.exit(options(saved))
  native = isTRUE(dbreg:::detect_backend(conn)[["is_databricks"]])
  table = basename(tempfile("ap_moment_test_"))
  big = paste0(table, "_big")
  drop = function(t) if (native) DBI::dbExecute(conn, paste("DROP TEMP TABLE IF EXISTS", t)) else DBI::dbRemoveTable(conn, t)
  on.exit(for (t in c(big, table)) try(drop(t), silent = TRUE), add = TRUE)
  # ODBC dbListTables omits session tables on native Databricks. SHOW TABLES
  # includes them, so cleanup assertions inspect the actual owned resources.
  tables = function() sort(if (native) DBI::dbGetQuery(conn, "SHOW TABLES")[["tableName"]] else DBI::dbListTables(conn))
  before = tables()
  set.seed(4082)
  n = 1200L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample.int(13, n, TRUE),
    h = sample.int(7, n, TRUE), j = sample.int(4, n, TRUE), k = sample(0:1, n, TRUE), all = 1,
    a = runif(n, .1, 5), cluster = sample.int(23, n, TRUE))
  f = sample.int(19, n, TRUE)
  d[["f"]] = paste0(c("A'", "B/", "Łódź_")[f %% 3 + 1], f)
  d[["endog"]] = .8 * d[["x"]] + .2 * d[["z"]] + rnorm(n)
  d[["x2"]] = d[["x"]]^2
  d[["y"]] = .4 * d[["endog"]] - .2 * d[["z"]] + .12 * d[["x"]]^2 + sin(f) + cos(d[["g"]]) + rnorm(n)
  if (native) DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", dbreg:::sql_values_frame(conn, d))) else
    DBI::dbWriteTable(conn, table, d, temporary = TRUE)
  # IDs that collide if converted to R doubles must remain distinct SQL keys.
  DBI::dbExecute(conn, dbreg:::backend_sql(conn, paste0("CREATE TEMP TABLE ", big,
    " AS SELECT * EXCLUDE(g),CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS g FROM ", table)))

  compare = function(fml, covariance, iv = FALSE, source = table) {
    if (verbose) message("[AP moments] ", paste(deparse(fml), collapse = ""), " ", deparse(covariance))
    fit = (if (iv) dbivreg else dbreg)(fml, conn = conn, table = source,
      weights = "a", vcov = covariance, strategy = "demean")
    ref = fixest::feols(fml, d, weights = ~a,
      vcov = if (identical(covariance, "hc1")) "hetero" else covariance,
      fixef.rm = "none", fixef.tol = 1e-10, ssc = fixest::ssc(K.fixef = "full"))
    ref_names = sub("^fit_", "", names(coef(ref)))
    order = match(names(coef(fit)), ref_names)
    expect_true(!anyNA(order))
    expect_true(all(abs(coef(fit) - coef(ref)[order]) <= 1e-7 + 1e-6 * abs(coef(ref)[order])))
    expect_true(all(abs(vcov(fit) - vcov(ref)[order, order]) <= 1e-8 + 1e-5 * abs(vcov(ref)[order, order])))
    expect_equal(nobs(fit), nrow(d))
    if (!iv) {
      expect_identical(fit[["ap_method"]], "bounded_moments")
      expect_true(fit[["ap_max_abs"]] <= 1e-10)
    }
    expect_identical(tables(), sort(c(before, table, big)))
    fit
  }
  for (fe in c("f + g", "f + g + h", "f + g + h + j + k + all")) {
    compare(as.formula(paste("y ~ x + z |", fe)), ~cluster)
  }
  fml = y ~ x + z + x2 | f + g + h
  for (covariance in list("iid", "hc1", ~g)) compare(fml, covariance)
  compare(fml, ~cluster, source = big)
  compare(y ~ z | f + g + h | endog ~ x, ~cluster, iv = TRUE)

  old = options(dbreg.ap_solver = "auto")
  tryCatch(compare(fml, ~cluster), finally = options(old))

  # All limits fall back to the established solver without changing estimates
  # or leaving temporary dictionaries behind. No wall-clock assertions in CI.
  reference = compare(y ~ x + z | f + g, ~cluster)
  for (cap in list(list(dbreg.duckdb_ap_max_edges = 1), list(dbreg.duckdb_ap_max_entries = 1),
                   list(dbreg.ap_moment_max_cells = 1), list(dbreg.glm_fe_dense_max = 0))) {
    old = options(cap)
    fit = tryCatch(dbreg(y ~ x + z | f + g, conn = conn, table = table,
      weights = "a", vcov = ~cluster, strategy = "demean"), finally = options(old))
    expect_true(fit[["ap_method"]] %in% c("row_tables", "window"))
    expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
    expect_true(max(abs(vcov(fit) - vcov(reference))) < 1e-8)
    expect_identical(tables(), sort(c(before, table, big)))
  }
  for (invalid in list(0, -1, 1.5, NA_real_, Inf)) {
    old = options(dbreg.ap_moment_max_cells = invalid)
    tryCatch(expect_error(dbreg(fml, conn = conn, table = table, strategy = "demean"),
      "positive whole numbers"), finally = options(old))
  }
  old = options(dbreg.ap_solver = "invalid")
  tryCatch(expect_error(dbreg(fml, conn = conn, table = table, strategy = "demean"), "ap_solver"), finally = options(old))

  # An unattainable accuracy request must fail, with all solve/fallback tables
  # released, rather than accepting an uncertified projected relation.
  old = options(dbreg.ap_tol = 1e-30, dbreg.ap_max_iter = 3L)
  tryCatch(expect_error(dbreg(fml, conn = conn, table = table, weights = "a", strategy = "demean"),
    "did not converge within 3 iterations"), finally = options(old))
  expect_identical(tables(), sort(c(before, table, big)))

  # Certification rejects lost rows and NULLs, even if all reported means
  # would otherwise be zero. It also detects a remaining group mean.
  certificate = function(value, expected = 1L) dbreg:::dbreg_projection_certificate(conn,
    paste0("(SELECT 1 AS f,1 AS g,1.0 AS w,", value, " AS x)"), c("f", "g"), "x", expected)
  expect_equal(certificate("0.0"), 0)
  expect_equal(certificate("0.0", 2L), Inf)
  expect_equal(certificate("CAST(NULL AS DOUBLE)"), Inf)
  expect_equal(certificate("0.25"), .25)
}
