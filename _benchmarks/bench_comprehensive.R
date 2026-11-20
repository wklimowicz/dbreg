#!/usr/bin/env Rscript
# 14-case comprehensive DuckDB (and optional Databricks) parity & performance suite
# Covers OLS, 1-3 way FE, composite FE, varying slopes, compressed OLS, IV 2SLS,
# Logit/Poisson dbglm, and dbbinsreg.
#
# Usage:
#   Rscript _benchmarks/bench_comprehensive.R
#
# Optional environment variables:
#   DBREG_BENCH_THREADS     DuckDB threads (default: 4)
#   DATABRICKS_SQL_PATH     Databricks SQL Warehouse HTTP path (optional; skipped if unset)
#   DBREG_BENCH_CSV         Output CSV path (default: _benchmarks/bench_comprehensive.csv)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dbreg)
})

THREADS  <- as.integer(Sys.getenv("DBREG_BENCH_THREADS", "4"))
DBX_PATH <- Sys.getenv("DATABRICKS_SQL_PATH", "")
CSV_OUT  <- Sys.getenv("DBREG_BENCH_CSV", "_benchmarks/bench_comprehensive.csv")

bench_uniform_sql <- function(id, seed) {
  key <- paste0("(CAST(", id, " AS BIGINT)+", as.integer(seed) * 104729L, ")")
  paste0("(CAST(MOD(MOD(", key, "*", key,
         ",2147483647)*48271,2147483647) AS DOUBLE)/2147483647.0)")
}

panel_sql <- function(n = 500000L) {
  x <- vapply(1:5, function(k) paste0(bench_uniform_sql("id", k), "-0.5 AS x", k), character(1))
  base <- paste0(
    "SELECT id, 1+MOD(id,10000) AS teacher,",
    " 1+MOD(MOD(id,10000)+CAST(FLOOR(id/10000) AS BIGINT),100) AS school,",
    " 1+MOD(CAST(FLOOR(id/10000) AS BIGINT),10) AS year,",
    " 0.5+", bench_uniform_sql("id", 25L), " AS w, ",
    paste(x, collapse = ", "),
    " FROM range(1,", n + 1L, ") t(id)"
  )
  eta <- "(-1.2+0.3*x1-0.2*x2+0.1*x3+0.2*sin(teacher*0.01)+0.1*cos(school*0.1))"
  paste0(
    "SELECT *, ",
    "0.2*x1-0.1*x2+0.08*x3+0.2*sin(teacher*0.01)+0.15*cos(school*0.1)+0.01*year+",
    bench_uniform_sql("id", 26L), "-0.5 AS yg, ",
    "0.5*x4 + 0.3*x2 + 0.2*(", bench_uniform_sql("id", 31L), "-0.5) AS d_endog, ",
    "CASE WHEN ", bench_uniform_sql("id", 27L), " < 1.0-EXP(-EXP(", eta, ")) THEN 1 ELSE 0 END AS y, ",
    "CAST(FLOOR(EXP(0.5 + 0.3*x1 - 0.2*x2 + 0.1*cos(school*0.1)) * (0.5 + ", bench_uniform_sql("id", 32L), ")) AS DOUBLE) AS y_cnt ",
    "FROM (", base, ") panel"
  )
}

discrete_sql <- function(n = 3100000L) {
  base <- paste0(
    "SELECT id, 1+MOD(id,296) AS district,",
    " 1+MOD(CAST(FLOOR(id/296) AS BIGINT),7) AS YEAR,",
    " MOD(CAST(FLOOR(id/2072) AS BIGINT),2) AS FTPTWK,",
    " MOD(CAST(FLOOR(id/4144) AS BIGINT),2) AS SEX,",
    " 20+MOD(CAST(FLOOR(id/8288) AS BIGINT),20) AS AGE",
    " FROM range(1,", n + 1L, ") t(id)"
  )
  paste0(
    "SELECT *, CASE WHEN ", bench_uniform_sql("id", 28L),
    " < 1.0/(1.0+EXP(-(-1.0+0.2*FTPTWK-0.15*SEX+0.01*(AGE-30)+0.1*sin(district*0.1)+0.02*YEAR)))",
    " THEN 1 ELSE 0 END AS permanent_wastage FROM (", base, ") discrete"
  )
}

cat(sprintf("Setting up DuckDB (%d threads)...\n", THREADS))
con_duck <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(dbDisconnect(con_duck, shutdown = TRUE), add = TRUE)
dbExecute(con_duck, sprintf("PRAGMA threads=%d", THREADS))
t0 <- proc.time()
dbExecute(con_duck, paste("CREATE TEMP TABLE bench_panel AS", panel_sql(500000L)))
dbExecute(con_duck, paste("CREATE TEMP TABLE bench_discrete AS", discrete_sql(3100000L)))
cat(sprintf("DuckDB materialized 500k panel + 3.1M discrete in %.2fs\n", (proc.time() - t0)[3]))

con_dbx <- NULL
if (nzchar(DBX_PATH) && requireNamespace("odbc", quietly = TRUE)) {
  con_dbx <- tryCatch(
    dbConnect(odbc::databricks(), httpPath = DBX_PATH, useNativeQuery = TRUE),
    error = function(e) { cat("Databricks connect failed:", conditionMessage(e), "\n"); NULL }
  )
}
if (!is.null(con_dbx)) {
  on.exit(dbDisconnect(con_dbx), add = TRUE)
  dbExecute(con_dbx, "SET use_cached_result = false")
  t0 <- proc.time()
  dbExecute(con_dbx, paste("CREATE TEMP TABLE bench_panel AS", panel_sql(500000L)))
  dbExecute(con_dbx, paste("CREATE TEMP TABLE bench_discrete AS", discrete_sql(3100000L)))
  cat(sprintf("Databricks materialized 500k panel + 3.1M discrete in %.2fs\n\n", (proc.time() - t0)[3]))
}

run_case <- function(name, fn) {
  cat(sprintf("--- Running: %s ---\n", name))
  run_one <- function(con, engine) {
    t0 <- proc.time()
    fit <- tryCatch(fn(con), error = function(e) structure(list(error = conditionMessage(e)), class = "bench_err"))
    elapsed <- (proc.time() - t0)[3]
    if (inherits(fit, "bench_err")) {
      cat(sprintf("  [%s] ERROR (%.2fs): %s\n", engine, elapsed, fit$error))
      return(list(elapsed = NA_real_, fit = NULL, error = fit$error, strategy = "ERROR"))
    }
    strat <- if (!is.null(fit$strategy)) fit$strategy else if (!is.null(fit$method)) fit$method else class(fit)[1]
    cat(sprintf("  [%s] %.3fs (strategy=%s)\n", engine, elapsed, strat))
    list(elapsed = elapsed, fit = fit, error = NULL, strategy = strat)
  }
  r_duck <- run_one(con_duck, "DuckDB")
  r_dbx  <- if (!is.null(con_dbx)) run_one(con_dbx, "Databricks") else list(elapsed = NA_real_, fit = NULL, strategy = "")

  coef_diff <- NA_real_
  se_diff <- NA_real_
  nobs_match <- NA
  if (!is.null(r_duck$fit) && !is.null(r_dbx$fit)) {
    if (inherits(r_duck$fit, "dbbinsreg")) {
      p_d <- r_duck$fit$points
      p_x <- r_dbx$fit$points
      if (!is.null(p_d) && !is.null(p_x) && nrow(p_d) == nrow(p_x)) {
        coef_diff <- max(abs(p_d$fit - p_x$fit), na.rm = TRUE)
        if ("se" %in% names(p_d) && "se" %in% names(p_x)) se_diff <- max(abs(p_d$se - p_x$se), na.rm = TRUE)
      }
      nobs_match <- identical(as.numeric(r_duck$fit$model$nobs_orig), as.numeric(r_dbx$fit$model$nobs_orig))
    } else {
      c_d <- coef(r_duck$fit); c_x <- coef(r_dbx$fit)
      common <- intersect(names(c_d), names(c_x))
      if (length(common)) coef_diff <- max(abs(c_d[common] - c_x[common]), na.rm = TRUE)
      se_d <- tryCatch(sqrt(diag(vcov(r_duck$fit))), error = function(e) NULL)
      se_x <- tryCatch(sqrt(diag(vcov(r_dbx$fit))), error = function(e) NULL)
      if (!is.null(se_d) && !is.null(se_x) && length(common)) se_diff <- max(abs(se_d[common] - se_x[common]), na.rm = TRUE)
      nobs_match <- identical(as.numeric(nobs(r_duck$fit)), as.numeric(nobs(r_dbx$fit)))
    }
  }
  data.frame(
    case = name,
    duckdb_s = round(r_duck$elapsed, 3),
    dbx_s = round(r_dbx$elapsed, 3),
    ratio_dbx_duck = round(r_dbx$elapsed / r_duck$elapsed, 2),
    duck_strategy = r_duck$strategy,
    dbx_strategy = r_dbx$strategy,
    max_coef_diff = coef_diff,
    max_se_diff = se_diff,
    nobs_match = nobs_match,
    stringsAsFactors = FALSE
  )
}

cases <- list(
  "01_OLS_HC1_500k" = function(con) dbreg(yg ~ x1 + x2 + x3, conn = con, table = "bench_panel", vcov = "hc1"),
  "02_1FE_10k_Cluster_500k" = function(con) dbreg(yg ~ x1 + x2 + x3 | teacher, conn = con, table = "bench_panel", vcov = ~school),
  "03_2FE_10kx100_Auto_500k" = function(con) dbreg(yg ~ x1 + x2 + x3 | teacher + school, conn = con, table = "bench_panel", vcov = "hc1"),
  "04_2FE_10kx100_IterAP_500k" = function(con) {
    op <- options(dbreg.ap_solver = "iterative"); on.exit(options(op))
    dbreg(yg ~ x1 + x2 + x3 | teacher + school, conn = con, table = "bench_panel", vcov = "hc1")
  },
  "05_3FE_10kx100x10_Auto_500k" = function(con) dbreg(yg ~ x1 + x2 + x3 | teacher + school + year, conn = con, table = "bench_panel", vcov = "hc1"),
  "06_CompositeFE_schoolXyear_500k" = function(con) dbreg(yg ~ x1 + x2 + x3 | school^year + teacher, conn = con, table = "bench_panel", vcov = "hc1"),
  "07_VaryingSlope_school_x3_500k" = function(con) dbreg(yg ~ x1 + x2 | school[x3], conn = con, table = "bench_panel", vcov = "hc1"),
  "08_Compress_OLS_3.1M" = function(con) dbreg(permanent_wastage ~ FTPTWK + SEX + AGE | district + YEAR, conn = con, table = "bench_discrete"),
  "09_IV_2SLS_2FE_500k" = function(con) dbivreg(yg ~ x1 + x2 | teacher + school | d_endog ~ x4, conn = con, table = "bench_panel", vcov = "hc1"),
  "10_GLM_Logit_Compress_3.1M" = function(con) dbglm(permanent_wastage ~ FTPTWK + SEX + AGE | district + YEAR, family = binomial(), conn = con, table = "bench_discrete"),
  "11_GLM_Logit_1FE_100_500k" = function(con) dbglm(y ~ x1 + x2 + x3 | school, family = binomial(), conn = con, table = "bench_panel"),
  "12_GLM_Poisson_1FE_100_500k" = function(con) dbglm(y_cnt ~ x1 + x2 + x3 | school, family = poisson(), conn = con, table = "bench_panel"),
  "13_Binscatter_NoFE_ContCtrl_500k" = function(con) dbbinsreg(yg ~ x1 + x2, conn = con, table = "bench_panel", nbins = 20, plot = FALSE),
  "14_Binscatter_1FE_DiscCtrl_500k" = function(con) dbbinsreg(yg ~ x1 + year | school, conn = con, table = "bench_panel", nbins = 20, plot = FALSE)
)

results <- do.call(rbind, Map(run_case, names(cases), cases))
cat("\n=== COMPREHENSIVE BENCHMARK SUMMARY ===\n")
print(results, row.names = FALSE)
write.csv(results, CSV_OUT, row.names = FALSE)
