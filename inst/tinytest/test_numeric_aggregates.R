# BIGINT aggregates must be decoded before entering ordinary double matrices.
# Exercise real DBI output rather than manufacturing an integer64 class.
local({
  modes = c("numeric", if (requireNamespace("bit64", quietly = TRUE)) "integer64")
  for (mode in modes) {
    conn = DBI::dbConnect(duckdb::duckdb(bigint = mode, config = list(threads = "2")))
    tryCatch({
      row = DBI::dbGetQuery(conn, paste(
        "SELECT CAST(10 AS BIGINT) AS sum_w, CAST(15 AS BIGINT) AS sum_w1,",
        "CAST(35 AS BIGINT) AS sum_w1_1, CAST(40 AS BIGINT) AS sum_w1_y,",
        "CAST(20 AS BIGINT) AS sum_wy, CAST(60 AS BIGINT) AS sum_wy_sq"))
      S = dbreg:::cross_moments(row, "x", "y", intercept = TRUE)
      expected = matrix(c(10, 15, 20, 15, 35, 40, 20, 40, 60), 3)
      expect_equal(as.numeric(S), as.numeric(expected))
      expect_equal(as.numeric(dbreg:::fit_ols(S, "y", c("(Intercept)", "x"))[["coef"]]), c(.8, .8))
      raw = DBI::dbGetQuery(conn, paste("SELECT CAST(3000000000 AS BIGINT) AS n,",
        "CAST(3000000001 AS BIGINT) AS dbreg_cell, CAST('9007199254740993' AS BIGINT) AS group_key"))
      decoded = dbreg:::decode_compressed_moments(raw)
      expect_identical(decoded[["n"]], 3e9)
      expect_identical(decoded[["dbreg_cell"]], 3000000001)
      expect_identical(decoded[["group_key"]], raw[["group_key"]])
    }, finally = DBI::dbDisconnect(conn, shutdown = TRUE))
  }
})

# Counts, weighted moments and compressed cell statistics feed several public
# estimators. All of them must agree with the same independent references.
local({
  saved = options(dbreg.ap_tol = 1e-10, dbreg.ap_max_iter = 200L)
  on.exit(options(saved), add = TRUE)
  fixest::setFixest_nthreads(2)
  set.seed(8103)
  n = 180L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(1:8, n, TRUE),
    h = sample(1:5, n, TRUE), cl = rep(1:12, each = 15), w = runif(n, .4, 2))
  d[["endog"]] = .8*d[["z"]] + .3*d[["x"]] + rnorm(n)
  d[["y"]] = .7 + .4*d[["x"]] + .8*d[["endog"]] + sin(d[["g"]]) + rnorm(n)
  modes = c("numeric", if (requireNamespace("bit64", quietly = TRUE)) "integer64")
  for (mode in modes) {
    conn = DBI::dbConnect(duckdb::duckdb(bigint = mode, config = list(threads = "2")))
    tryCatch({
      DBI::dbWriteTable(conn, "aggregate_source", d)
      # FE labels above 2^53 are exact keys, not numeric moment columns.
      DBI::dbExecute(conn, paste(
        "CREATE TEMP TABLE aggregate_keys AS SELECT *,",
        "CAST('9007199254740992' AS BIGINT) + g AS large_g FROM aggregate_source"))
      if (mode == "numeric") DBI::dbExecute(conn, paste(
        "CREATE OR REPLACE TEMP TABLE aggregate_keys AS SELECT * EXCLUDE (large_g),",
        "CAST(large_g AS VARCHAR) AS large_g FROM aggregate_keys"))
      for (strategy in c("moments", "demean", "compress", "mundlak")) {
        has_fe = strategy %in% c("demean", "mundlak")
        fml = if (has_fe) y ~ x | g else y ~ x
        for (weighted in c(FALSE, TRUE)) {
          if (strategy == "mundlak" && weighted) next # weighted CRE is a separate contract
          for (kind in list("iid", "hc1", ~cl)) {
            args = list(fml, conn = conn, table = "aggregate_source", strategy = strategy,
              vcov = kind, weights = if (weighted) "w" else NULL)
            fit = tryCatch(do.call(dbreg, args), error = function(e) e)
            expect_false(inherits(fit, "error"), info = paste(mode, strategy, weighted, deparse(kind)))
            if (inherits(fit, "error")) next
            if (strategy == "mundlak") {
              augmented = d
              augmented[["group_mean"]] = ave(d[["x"]], d[["g"]])
              ref = lm(y ~ x + group_mean, augmented)
              X = model.matrix(ref); bread = solve(crossprod(X)); scores = X*residuals(ref)
              V = if (is.character(kind) && kind == "iid") vcov(ref) else {
                meat = if (is.character(kind)) crossprod(scores) else crossprod(rowsum(scores, d[["cl"]]))
                correction = if (is.character(kind)) n/(n-ncol(X)) else 12/11*(n-1)/(n-ncol(X))
                correction*bread%*%meat%*%bread
              }
              dimnames(V) = list(colnames(X), colnames(X))
              expect_equal(as.numeric(coef(fit)["x"]), as.numeric(coef(ref)["x"]), tolerance = 1e-8)
              expect_equal(as.numeric(vcov(fit)["x", "x"]), as.numeric(V["x", "x"]), tolerance = 1e-8)
            } else {
              ref = fixest::feols(fml, d, weights = if (weighted) ~w else NULL, vcov = kind,
                fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
              terms = names(coef(ref))
              expect_equal(as.numeric(coef(fit)[terms]), as.numeric(coef(ref)), tolerance = 1e-8)
              expect_equal(as.numeric(as.matrix(vcov(fit))[terms, terms, drop = FALSE]), as.numeric(vcov(ref)), tolerance = 1e-8)
            }
            expect_identical(fit[["nobs_orig"]], as.numeric(n))
            expect_true(is.numeric(fit[["df_residual"]]) && !inherits(fit[["df_residual"]], "integer64"))
          }
        }
      }
      for (strategy in c("moments", "compress", "demean")) {
        fml = if (strategy == "demean") y ~ x | g + h | endog ~ z else y ~ x | endog ~ z
        for (kind in list("iid", "hc1", ~cl)) {
          fit = tryCatch(dbivreg(fml, conn = conn, table = "aggregate_source", strategy = strategy,
            vcov = kind, weights = "w"), error = function(e) e)
          expect_false(inherits(fit, "error"), info = paste(mode, strategy, "IV", deparse(kind)))
          if (inherits(fit, "error")) next
          ref = fixest::feols(fml, d, weights = ~w, vcov = kind, fixef.rm = "none",
            ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
          b = coef(ref); names(b) = sub("^fit_", "", names(b))
          V = vcov(ref); dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
          terms = names(coef(fit))
          expect_equal(as.numeric(coef(fit)), as.numeric(b[terms]), tolerance = 1e-8)
          expect_true(max(abs(as.numeric(vcov(fit)) - as.numeric(V[terms, terms]))) < 1e-8,
            info = paste(mode, strategy, "IV covariance", deparse(kind)))
          expect_identical(fit[["nobs_orig"]], as.numeric(n))
        }
      }
      # Disabling missing-value filtering selects the cluster-by-cell
      # fallback rather than the SQL score cross-product path. No values here
      # are missing, so both paths have the same independent clustered oracle.
      for (iv in c(FALSE, TRUE)) {
        fml = if (iv) y ~ x | endog ~ z else y ~ x
        estimator = if (iv) dbivreg else dbreg
        fit = estimator(fml, conn = conn, table = "aggregate_source", strategy = "compress",
          vcov = ~cl, drop_missings = FALSE)
        ref = fixest::feols(fml, d, vcov = ~cl, fixef.rm = "none")
        b = coef(ref); names(b) = sub("^fit_", "", names(b))
        V = vcov(ref); dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
        terms = names(coef(fit))
        expect_equal(as.numeric(coef(fit)), as.numeric(b[terms]), tolerance = 1e-8)
        expect_true(max(abs(as.numeric(vcov(fit)) - as.numeric(V[terms, terms]))) < 1e-8)
        expect_identical(fit[["nobs_orig"]], as.numeric(n))
      }
      # Converting aggregates must not round the original cell keys.
      ref = fixest::feols(y ~ x | g, d, vcov = ~cl, fixef.rm = "none",
        ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
      fit = tryCatch(dbreg(y ~ x | large_g, conn = conn, table = "aggregate_keys",
        strategy = "compress", vcov = ~cl), error = function(e) e)
      expect_false(inherits(fit, "error"))
      if (!inherits(fit, "error")) {
        expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
        V = as.matrix(fit[["vcov"]]); dimnames(V) = rep(list(rownames(fit[["coeftable"]])), 2L)
        expect_equal(as.numeric(V["x", "x"]), as.numeric(vcov(ref)), tolerance = 1e-8)
        expect_identical(length(unique(fit[["data"]][["large_g"]])), 8L)
      }
      expect_identical(sort(DBI::dbListTables(conn)), c("aggregate_keys", "aggregate_source"))
      expect_true(DBI::dbIsValid(conn))
    }, finally = DBI::dbDisconnect(conn, shutdown = TRUE))
  }
})
