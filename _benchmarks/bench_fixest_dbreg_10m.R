#!/usr/bin/env Rscript
# 10M-row fixest + dbreg benchmark suite:
#   Compares fixest (in-memory R) vs. dbreg {fork, optional upstream} x {DuckDB, optional Databricks}
#
# Usage:
#   Rscript _benchmarks/bench_fixest_dbreg_10m.R
#
# Optional environment variables:
#   DBREG_BENCH_N           Row count (default: 10000000)
#   DBREG_BENCH_THREADS     Threads for fixest & DuckDB (default: 8)
#   DBREG_UPSTREAM_LIB      Library path containing upstream dbreg (optional; skipped if unset)
#   DATABRICKS_SQL_PATH     Databricks SQL Warehouse HTTP path (optional; skipped if unset)
#   DBREG_BENCH_CSV         Output CSV path (default: _benchmarks/bench_fixest_dbreg_10m.csv)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(fixest)
})

N        <- as.integer(Sys.getenv("DBREG_BENCH_N", "10000000"))
THREADS  <- as.integer(Sys.getenv("DBREG_BENCH_THREADS", "8"))
UP_LIB   <- Sys.getenv("DBREG_UPSTREAM_LIB", "")
DBX_PATH <- Sys.getenv("DATABRICKS_SQL_PATH", "")
CSV_OUT  <- Sys.getenv("DBREG_BENCH_CSV", "_benchmarks/bench_fixest_dbreg_10m.csv")

setFixest_nthreads(THREADS)

# Deterministic SQL hash uniform in [0, 1) that evaluates identically on DuckDB & Databricks
bench_u <- function(id, seed) {
  key <- paste0("(CAST(", id, " AS BIGINT)+", as.integer(seed) * 104729L, ")")
  paste0("(CAST(MOD(MOD(", key, "*", key,
         ",2147483647)*48271,2147483647) AS DOUBLE)/2147483647.0)")
}

# 1. fixest-style panel: Easy (independent FEs) + Difficult (weak-mobility / correlated FEs)
#    plus OLS, Poisson, and Logit outcomes
fixest_sql <- paste0(
  "SELECT ",
  "  id, ",
  "  1 + MOD(id, 10000) AS fe1_e, ",
  "  1 + MOD(CAST(FLOOR(id * 7 / 13) AS BIGINT), 1000) AS fe2_e, ",
  "  1 + MOD(CAST(FLOOR(id * 11 / 17) AS BIGINT), 100) AS fe3_e, ",
  "  1 + MOD(id, 10000) AS fe1_d, ",
  "  1 + MOD(CAST(FLOOR(MOD(id, 10000) / 10) AS BIGINT) + MOD(CAST(FLOOR(id / 10000) AS BIGINT), 3), 1000) AS fe2_d, ",
  "  1 + MOD(CAST(FLOOR(MOD(id, 10000) / 100) AS BIGINT) + MOD(CAST(FLOOR(id / 10000) AS BIGINT), 2), 100) AS fe3_d, ",
  "  1 + MOD(id, 500) AS g1, ",
  "  1 + MOD(CAST(FLOOR(id / 500) AS BIGINT), 50) AS g2, ",
  "  ", bench_u("id", 1L), " - 0.5 AS x1, ",
  "  ", bench_u("id", 2L), " - 0.5 AS x2, ",
  "  ", bench_u("id", 3L), " - 0.5 AS x3, ",
  "  0.5 * (", bench_u("id", 1L), " - 0.5) - 0.3 * (", bench_u("id", 2L), " - 0.5) + 0.2 * (", bench_u("id", 3L), " - 0.5) ",
  "    + 0.2 * SIN(MOD(id, 10000) * 0.01) + 0.15 * COS(MOD(CAST(FLOOR(id / 10000) AS BIGINT), 100) * 0.1) + (", bench_u("id", 4L), " - 0.5) AS y_ols, ",
  "  CAST(FLOOR(EXP(0.2 + 0.4 * (", bench_u("id", 1L), " - 0.5) - 0.2 * (", bench_u("id", 2L), " - 0.5) + 0.1 * SIN(MOD(id, 500) * 0.05) + 0.5 * ", bench_u("id", 5L), ")) AS DOUBLE) AS y_pois, ",
  "  CASE WHEN ", bench_u("id", 6L), " < 1.0 / (1.0 + EXP(-(-0.2 + 0.5 * (", bench_u("id", 1L), " - 0.5) - 0.3 * (", bench_u("id", 2L), " - 0.5) + 0.15 * SIN(MOD(id, 500) * 0.05)))) THEN 1.0 ELSE 0.0 END AS y_bin ",
  "FROM range(1, ", N + 1L, ") t(id)"
)

# 2. dbreg discrete 2-way panel (compress benchmark: 2,000 units x 5,000 periods, K1 = K2 = 4)
dbreg_disc_sql <- paste0(
  "SELECT ",
  "  1 + MOD(id - 1, 2000) AS unit, ",
  "  1 + CAST(FLOOR((id - 1) / 2000) AS BIGINT) AS period, ",
  "  1 + MOD(CAST(FLOOR(MOD(id, 997) * 4 / 997) AS BIGINT), 4) AS d1, ",
  "  1 + MOD(CAST(FLOOR(MOD(id * 3, 991) * 4 / 991) AS BIGINT), 4) AS d2, ",
  "  0.5 * (1 + MOD(CAST(FLOOR(MOD(id, 997) * 4 / 997) AS BIGINT), 4)) ",
  "    - 0.3 * (1 + MOD(CAST(FLOOR(MOD(id * 3, 991) * 4 / 991) AS BIGINT), 4)) ",
  "    + 0.2 * SIN(MOD(id - 1, 2000) * 0.01) + 0.1 * COS(CAST(FLOOR((id - 1) / 2000) AS BIGINT) * 0.01) ",
  "    + (", bench_u("id", 9L), " - 0.5) AS y ",
  "FROM range(1, ", N + 1L, ") t(id)"
)

cat(sprintf("=== Setting up %.1fM-row fixest + dbreg benchmark (%d threads) ===\n", N / 1e6, THREADS))
flush.console()

con_duck <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(dbDisconnect(con_duck, shutdown = TRUE), add = TRUE)
dbExecute(con_duck, sprintf("PRAGMA threads=%d", THREADS))

t0 <- proc.time()[["elapsed"]]
dbExecute(con_duck, paste("CREATE TABLE fixest_10m AS", fixest_sql))
dbExecute(con_duck, paste("CREATE TABLE dbreg_disc_10m AS", dbreg_disc_sql))
cat(sprintf("DuckDB materialized both tables in %.2fs\n", proc.time()[["elapsed"]] - t0))

t0 <- proc.time()[["elapsed"]]
df_fixest <- dbGetQuery(con_duck, "SELECT * FROM fixest_10m")
df_disc   <- dbGetQuery(con_duck, "SELECT * FROM dbreg_disc_10m")
cat(sprintf("Loaded into R data.frame for fixest in %.2fs\n", proc.time()[["elapsed"]] - t0))
flush.console()

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
  t0 <- proc.time()[["elapsed"]]
  dbExecute(con_dbx, paste("CREATE TEMP TABLE fixest_10m AS", fixest_sql))
  dbExecute(con_dbx, paste("CREATE TEMP TABLE dbreg_disc_10m AS", dbreg_disc_sql))
  cat(sprintf("Databricks materialized both TEMP tables in %.2fs\n\n", proc.time()[["elapsed"]] - t0))
  flush.console()
}

has_upstream <- nzchar(UP_LIB) && dir.exists(UP_LIB)
default_libs <- .libPaths()
use_version <- function(which_ver) {
  if ("package:dbreg" %in% search()) detach("package:dbreg", unload = TRUE, character.only = TRUE)
  if ("dbreg" %in% loadedNamespaces()) unloadNamespace("dbreg")
  if (which_ver == "fork") .libPaths(default_libs) else .libPaths(c(UP_LIB, default_libs))
  suppressPackageStartupMessages(library(dbreg))
}

cases <- list(
  list(
    name = "01_dbreg_Discrete_2FE_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y ~ d1 + d2 | unit + period, data = df_disc, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y ~ d1 + d2 | unit + period, conn = con, table = "dbreg_disc_10m")
  ),
  list(
    name = "02_fixest_OLS_1FE_Easy_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y_ols ~ x1 + x2 + x3 | fe1_e, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y_ols ~ x1 + x2 + x3 | fe1_e, conn = con, table = "fixest_10m")
  ),
  list(
    name = "03_fixest_OLS_2FE_Easy_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y_ols ~ x1 + x2 + x3 | fe1_e + fe2_e, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y_ols ~ x1 + x2 + x3 | fe1_e + fe2_e, conn = con, table = "fixest_10m")
  ),
  list(
    name = "04_fixest_OLS_3FE_Easy_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y_ols ~ x1 + x2 + x3 | fe1_e + fe2_e + fe3_e, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y_ols ~ x1 + x2 + x3 | fe1_e + fe2_e + fe3_e, conn = con, table = "fixest_10m")
  ),
  list(
    name = "05_fixest_OLS_2FE_Difficult_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y_ols ~ x1 + x2 + x3 | fe1_d + fe2_d, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y_ols ~ x1 + x2 + x3 | fe1_d + fe2_d, conn = con, table = "fixest_10m")
  ),
  list(
    name = "06_fixest_OLS_3FE_Difficult_10M",
    has_upstream = TRUE,
    fixest_fn = function() fixest::feols(y_ols ~ x1 + x2 + x3 | fe1_d + fe2_d + fe3_d, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbreg(y_ols ~ x1 + x2 + x3 | fe1_d + fe2_d + fe3_d, conn = con, table = "fixest_10m")
  ),
  list(
    name = "07_fixest_Poisson_2FE_10M",
    has_upstream = FALSE,
    fixest_fn = function() fixest::fepois(y_pois ~ x1 + x2 | g1 + g2, data = df_fixest, fixef.rm = "none"),
    dbreg_fn  = function(con) dbglm(y_pois ~ x1 + x2 | g1 + g2, family = poisson(), conn = con, table = "fixest_10m")
  ),
  list(
    name = "08_fixest_Logit_2FE_10M",
    has_upstream = FALSE,
    fixest_fn = function() fixest::feglm(y_bin ~ x1 + x2 | g1 + g2, data = df_fixest, family = binomial(), fixef.rm = "none"),
    dbreg_fn  = function(con) dbglm(y_bin ~ x1 + x2 | g1 + g2, family = binomial(), conn = con, table = "fixest_10m")
  )
)

run_timed <- function(label, fn) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(fn(), error = function(e) structure(list(error = conditionMessage(e)), class = "bench_err"))
  dt <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "bench_err")) {
    cat(sprintf("  [%-20s] ERROR (%.2fs): %s\n", label, dt, res$error))
    flush.console()
    return(list(sec = NA_real_, coef = NULL, err = res$error))
  }
  cf <- as.numeric(coef(res))
  cat(sprintf("  [%-20s] %7.3fs\n", label, dt))
  flush.console()
  list(sec = dt, coef = cf, err = "")
}

rows <- list()
for (ci in cases) {
  cat(sprintf("--- %s ---\n", ci$name))
  flush.console()
  r_fe <- run_timed("fixest (R)", ci$fixest_fn)
  use_version("fork")
  r_fd <- run_timed("fork | DuckDB", function() ci$dbreg_fn(con_duck))
  r_ud <- if (has_upstream && ci$has_upstream) { use_version("upstream"); run_timed("upst | DuckDB", function() ci$dbreg_fn(con_duck)) } else list(sec = NA_real_, coef = NULL, err = "")
  r_fb <- if (!is.null(con_dbx)) { use_version("fork"); run_timed("fork | Databricks", function() ci$dbreg_fn(con_dbx)) } else list(sec = NA_real_, coef = NULL, err = "")
  r_ub <- if (!is.null(con_dbx) && has_upstream && ci$has_upstream) { use_version("upstream"); run_timed("upst | Databricks", function() ci$dbreg_fn(con_dbx)) } else list(sec = NA_real_, coef = NULL, err = "")

  d_fe <- if (!is.null(r_fe$coef) && !is.null(r_fd$coef)) max(abs(r_fe$coef - r_fd$coef)) else NA_real_
  d_db <- if (!is.null(r_fd$coef) && !is.null(r_fb$coef)) max(abs(r_fd$coef - r_fb$coef)) else NA_real_

  rows[[length(rows) + 1L]] <- data.frame(
    case = ci$name,
    fixest_s = round(r_fe$sec, 3),
    fork_duck_s = round(r_fd$sec, 3),
    upst_duck_s = round(r_ud$sec, 3),
    fork_dbx_s  = round(r_fb$sec, 3),
    upst_dbx_s  = round(r_ub$sec, 3),
    diff_vs_fixest = signif(d_fe, 3),
    diff_duck_dbx  = signif(d_db, 3),
    stringsAsFactors = FALSE
  )
  write.csv(do.call(rbind, rows), CSV_OUT, row.names = FALSE)
}

cat("\n=== DONE ===\n")
print(do.call(rbind, rows), row.names = FALSE)
