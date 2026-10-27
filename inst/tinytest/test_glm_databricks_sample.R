library(DBI)
source(system.file("tinytest/helper_glm_reference.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_glm_databricks_pruning.R", package = "dbreg"), local = TRUE)
local({
  live = identical(tolower(Sys.getenv("DBREG_TEST_DATABRICKS")), "true")
  if (live) {
    native = dbConnect(odbc::databricks(),
      driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
      httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
    conn = native
    on.exit(dbDisconnect(native))
    dbExecute(conn, "SET use_cached_result=false")
  } else {
    native = dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2")))
    conn = dbglm_test_databricks_connection(native)
    on.exit(dbDisconnect(native, shutdown = TRUE))
  }
  dbglm_test_databricks_pruning(conn)
  source = basename(tempfile("glm_sample_"))
  input = paste0(source, "_input")
  retained = paste0(source, "_retained")
  failed = paste0(source, "_failed")
  clean_table = paste0(source, "_clean")
  on.exit(for (table in c(clean_table, retained, input, source)) try(dbExecute(conn, paste("DROP TEMP TABLE", table)), silent = TRUE), add = TRUE, after = FALSE)
  set.seed(584)
  d = data.frame(id = seq_len(480), x = 10.8 + rnorm(480), f = rep(1:4, 120),
    g = rep(1:6, each = 80), a = 1, o = 0, cl = rep(1:12, each = 40))
  d[["y"]] = rbinom(480, 1, plogis(-.5 + .2 * (d[["x"]] - 10.8)))
  extra = d[1:4, ]
  extra[["id"]] = 481:484; extra[["f"]] = c(5, 5, 6, 6)
  extra[["g"]] = c(1, 7, 7, 8); extra[["y"]] = c(0, 1, 0, 1)
  d = rbind(d, extra)
  d[["y"]][d[["f"]] == 3] = 0L # Leaves a hole in the active FE ranks.
  for (j in 1:4) d[[paste0("h", j)]] = sample(1:3, nrow(d), TRUE)
  dbExecute(conn, paste("CREATE TEMP TABLE", source, "AS", dbreg:::sql_values_frame(conn, d)))
  dbExecute(conn, paste0("CREATE TEMP TABLE ", input,
    " AS SELECT ", dbreg:::sql_select_except(conn, c("x", "o", "a", "cl")), ",CASE WHEN id=1 THEN NULL ELSE x END AS x,",
    "CASE WHEN id=2 THEN CAST('Infinity' AS DOUBLE) ELSE o END AS o,",
    "CASE WHEN id=3 THEN 0 ELSE a END AS a,CASE WHEN id=4 THEN NULL ELSE cl END AS cl FROM ", source))
  formula = y ~ x | f + g + h1 + h2 + h3 + h4
  eligible = d[-(1:4), ]
  ref = dbglm_reference_feglm(formula, eligible, family = binomial("cloglog"),
    weights = ~a, offset = ~o, vcov = ~cl, fixef.rm = "perfect_fit")
  fit = dbglm(formula, conn = conn, table = input, family = binomial("cloglog"),
    weights = "a", offset = "o", vcov = ~cl, row_id = "id", sample_table = retained)
  keys = sort(as.numeric(dbGetQuery(conn, paste("SELECT id FROM", retained))[["id"]]))
  expect_equal(keys, sort(as.numeric(eligible[["id"]][fixest::obs(ref)])))
  expect_equal(nobs(fit), nobs(ref))
  expect_equal(fit[["nobs_dropped"]], nrow(d) - length(keys))
  expect_equal(fit[["nobs_eligible"]], nrow(eligible))
  expect_true(fit[["pruning_rounds"]] > 1L)
  expect_true(fit[["score_max"]] <= fit[["control"]][["score_tol"]])
  expect_equal(coef(fit), coef(ref), tolerance = 1e-7)
  expect_true(max(abs(vcov(fit) - vcov(ref))) < 1e-7)
  expect_equal(fit[["n_params_vcov"]], attr(vcov(ref, attr = TRUE), "df.K"))
  expect_error(dbglm(family = stats::binomial("cloglog"), formula, conn = conn, table = input, weights = "a", offset = "o", vcov = ~cl,
    fixef_rm = "none"), "separation")
  expect_error(dbglm(family = stats::binomial("cloglog"), formula, conn = conn, table = input, weights = "a", offset = "o", vcov = ~cl,
    sample_table = failed, control = dbglm_control(glm_max_iter = 1)), "converge")
  expect_false(dbExistsTable(conn, failed))
  expect_equal(as.numeric(dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", input))[["n"]]), nrow(d))
  expect_true(dbIsValid(conn))
  # Exercise the skipped-pruning path end to end with a compound FE, against
  # fixest's independently constructed sample and full clustered covariance.
  clean = d[d[["id"]] <= 480 & d[["f"]] != 3, ]
  dbExecute(conn, paste("CREATE TEMP TABLE", clean_table, "AS SELECT * FROM", source, "WHERE id<=480 AND f!=3"))
  expect_false(dbreg:::dbglm_pruning_needed(conn, clean_table,
    list(c("f", "h1"), "g", "h2", "h3", "h4"), "y", "perfect_fit", "binomial"))
  composite = y ~ x | f^h1 + g + h2 + h3 + h4
  ref_clean = dbglm_reference_feglm(composite, clean, family = binomial("cloglog"), vcov = ~g, fixef.rm = "none")
  fit_clean = dbglm(composite, conn = conn, table = clean_table, family = binomial("cloglog"), vcov = ~g)
  expect_equal(fit_clean[["pruning_rounds"]], 0L)
  expect_equal(nobs(fit_clean), nobs(ref_clean))
  expect_equal(coef(fit_clean), coef(ref_clean), tolerance = 1e-7)
  expect_true(max(abs(vcov(fit_clean) - vcov(ref_clean))) < 1e-7)
  expect_equal(fit_clean[["n_params_vcov"]], attr(vcov(ref_clean, attr = TRUE), "df.K"))
})
