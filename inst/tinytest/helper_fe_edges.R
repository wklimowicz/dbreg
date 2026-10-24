# Small, deterministic stress cases shared by DuckDB and native Databricks.
# The reference always uses the original data, independently of dbreg's sample.
dbreg_test_fe_edge_data = function() {
  set.seed(20261024)
  make_data = function(n) {
    d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
      f = sample.int(12, n, TRUE), g = sample.int(9, n, TRUE),
      a = 1, o = 0, cl = sample.int(17, n, TRUE))
    d
  }
  disconnected = make_data(1800)
  component = rep(0:2, each = 600)
  disconnected[["f"]] = component * 4 + sample.int(4, 1800, TRUE)
  disconnected[["g"]] = component * 3 + sample.int(3, 1800, TRUE)

  weighted = make_data(1800)
  weighted[["a"]] = sample(c(1, 2, 10, 100), 1800, TRUE, c(.6, .2, .15, .05))
  weighted[["o"]] = runif(1800, -1.5, .8)

  redundant = make_data(1800)
  redundant[["duplicate"]] = 2 * redundant[["x"]]
  redundant[["absorbed"]] = sin(redundant[["f"]])

  rare = make_data(6000)
  rare[["f"]] = sample.int(12, 6000, TRUE, prob = seq_len(12))
  rare[["g"]] = (rare[["f"]] - 1L) %% 9L + 1L
  move = seq_len(900)
  rare[["g"]][move] = sample.int(9, length(move), TRUE)
  rare[["o"]] = runif(6000, -.6, .6)

  cases = list(disconnected = disconnected, unequal_weights = weighted,
    redundant_columns = redundant, rare_correlated = rare)
  lapply(names(cases), function(name) {
    d = cases[[name]]
    linear = .3 * d[["x"]] - .2 * d[["z"]] + .15 * sin(d[["f"]]) + .1 * cos(d[["g"]])
    eta = (if (name == "rare_correlated") -4.2 else -1.2) + linear + d[["o"]]
    d[["binary"]] = rbinom(nrow(d), 1, binomial("cloglog")[["linkinv"]](eta))
    d[["continuous"]] = linear + rnorm(nrow(d))
    rhs = if (name == "redundant_columns") "x + duplicate + z + absorbed" else "x + z"
    list(name = name, data = d, glm = as.formula(paste("binary ~", rhs, "| f + g")),
      ols = as.formula(paste("continuous ~", rhs, "| f + g")))
  })
}

dbreg_test_fe_edges = function(conn, estimator = c("glm", "ols"), verbose = FALSE) {
  estimator = match.arg(estimator)
  is_databricks = isTRUE(dbreg:::detect_backend(conn)[["is_databricks"]])
  saved = options(dbreg.ap_tol = 1e-10, dbreg.ap_max_iter = 1000L)
  on.exit(options(saved))
  for (case in dbreg_test_fe_edge_data()) {
    d = case[["data"]]
    name = paste(estimator, case[["name"]])
    if (verbose) message("[FE edge test] ", name)
    table = basename(tempfile("fe_edge_"))
    if (is_databricks) {
      DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", dbreg:::dbglm_values_sql(conn, d)))
    } else DBI::dbWriteTable(conn, table, d, temporary = TRUE)
    tryCatch({
      # An iid and an externally clustered covariance exercise both the
      # final-information and score-aggregation passes on each difficult design.
      for (kind in c("iid", "cluster")) {
        covariance = if (kind == "cluster") ~cl else "iid"
        if (estimator == "glm") {
          ref = dbglm_reference_feglm(case[["glm"]], d, family = binomial("cloglog"),
            weights = ~a, offset = ~o, vcov = covariance, fixef.rm = "none")
          got = dbglm(case[["glm"]], conn = conn, table = table, family = binomial("cloglog"),
            weights = "a", offset = "o", vcov = covariance, fixef_rm = "none")
          expect_true(got[["converged"]] && got[["score_max"]] <= got[["control"]][["score_tol"]], info = name)
          expect_true(got[["projection_residual"]] <= got[["control"]][["fe_tol"]], info = name)
          expect_true(abs(as.numeric(logLik(got) - logLik(ref))) / nobs(ref) < 1e-9, info = name)
          expect_equal(got[["n_params_vcov"]], attr(vcov(ref, attr = TRUE), "df.K"), info = name)
        } else {
          ref = fixest::feols(case[["ols"]], d, weights = ~a, vcov = covariance,
            fixef.rm = "none", fixef.tol = 1e-10, ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
          got = dbreg(case[["ols"]], conn = conn, table = table,
            weights = "a", vcov = covariance, strategy = "demean")
        }
        expect_identical(names(coef(got)), names(coef(ref)), info = name)
        expect_true(all(is.finite(coef(got))) && all(is.finite(vcov(got))), info = name)
        expect_true(all(abs(coef(got) - coef(ref)) <= 1e-7 + 1e-6 * abs(coef(ref))), info = name)
        expect_true(all(abs(vcov(got) - vcov(ref)) <= 1e-8 + 1e-5 * abs(vcov(ref))), info = name)
        expect_equal(nobs(got), nobs(ref), info = name)
        if (case[["name"]] == "redundant_columns") {
          expect_identical(names(coef(got)), c("x", "z"), info = name)
        }
      }
      expect_equal(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]]), nrow(d), info = name)
    }, finally = if (is_databricks) DBI::dbExecute(conn, paste("DROP TEMP TABLE", table)) else DBI::dbRemoveTable(conn, table))
  }
}
