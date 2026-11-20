#!/usr/bin/env Rscript
# 12 GB (100,000,000 rows x 15 cols) out-of-core benchmark:
#   Disk-backed DuckDB (8GB memory cap) vs. optional Databricks, Fork vs. optional Upstream
#
# Usage:
#   Rscript _benchmarks/bench_12gb.R
#
# Optional environment variables:
#   DBREG_BENCH_N           Row count (default: 100000000)
#   DBREG_BENCH_THREADS     DuckDB threads (default: 4)
#   DBREG_BENCH_MEM_LIMIT   DuckDB memory_limit (default: 8GB)
#   DBREG_BENCH_DIR         Directory for on-disk .duckdb file & spill tmp (default: tempdir())
#   DBREG_UPSTREAM_LIB      Library path containing upstream dbreg (optional; skipped if unset)
#   DATABRICKS_SQL_PATH     Databricks SQL Warehouse HTTP path (optional; skipped if unset)
#   DBREG_BENCH_CSV         Output CSV path (default: _benchmarks/bench_12gb_results.csv)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

N_ROWS   <- as.numeric(Sys.getenv("DBREG_BENCH_N", "100000000"))
THREADS  <- as.integer(Sys.getenv("DBREG_BENCH_THREADS", "4"))
MEM_LIM  <- Sys.getenv("DBREG_BENCH_MEM_LIMIT", "8GB")
WORK_DIR <- Sys.getenv("DBREG_BENCH_DIR", tempdir())
UP_LIB   <- Sys.getenv("DBREG_UPSTREAM_LIB", "")
DBX_PATH <- Sys.getenv("DATABRICKS_SQL_PATH", "")
CSV_OUT  <- Sys.getenv("DBREG_BENCH_CSV", "_benchmarks/bench_12gb_results.csv")

duck_file <- file.path(WORK_DIR, "bench_12gb.duckdb")
duck_tmp  <- file.path(WORK_DIR, "duckdb_12gb_tmp")
unlink(c(duck_file, paste0(duck_file, ".wal"), duck_tmp), recursive = TRUE)
dir.create(duck_tmp, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(c(duck_file, paste0(duck_file, ".wal"), duck_tmp), recursive = TRUE), add = TRUE)

bench_uniform_sql <- function(id, seed) {
  key <- paste0("(CAST(", id, " AS BIGINT)+", as.integer(seed) * 104729L, ")")
  paste0("(CAST(MOD(MOD(", key, "*", key,
         ",2147483647)*48271,2147483647) AS DOUBLE)/2147483647.0)")
}

table_100m_sql <- function(n = 100000000) {
  x <- vapply(1:3, function(k) paste0(bench_uniform_sql("id", k), "-0.5 AS x", k), character(1))
  base <- paste0(
    "SELECT id,",
    " 1+MOD(id,10000) AS teacher,",
    " 1+MOD(MOD(id,10000)+CAST(FLOOR(id/10000) AS BIGINT),100) AS school,",
    " 1+MOD(CAST(FLOOR(id/10000) AS BIGINT),10) AS year,",
    " 1+MOD(id,296) AS district,",
    " 1+MOD(CAST(FLOOR(id/296) AS BIGINT),7) AS YEAR_D,",
    " MOD(CAST(FLOOR(id/2072) AS BIGINT),2) AS FTPTWK,",
    " MOD(CAST(FLOOR(id/4144) AS BIGINT),2) AS SEX,",
    " 20+MOD(CAST(FLOOR(id/8288) AS BIGINT),20) AS AGE,",
    " 0.5+", bench_uniform_sql("id", 25L), " AS w, ",
    paste(x, collapse = ", "),
    sprintf(" FROM range(1, %.0f) t(id)", n + 1)
  )
  paste0(
    "SELECT *, ",
    "0.2*x1-0.1*x2+0.08*x3+0.2*sin(teacher*0.01)+0.15*cos(school*0.1)+0.01*year+",
    bench_uniform_sql("id", 26L), "-0.5 AS yg, ",
    "CASE WHEN ", bench_uniform_sql("id", 28L),
    " < 1.0/(1.0+EXP(-(-1.0+0.2*FTPTWK-0.15*SEX+0.01*(AGE-30)+0.1*sin(district*0.1)+0.02*YEAR_D)))",
    " THEN 1.0 ELSE 0.0 END AS permanent_wastage FROM (", base, ") panel100m"
  )
}

cat(sprintf("=== Setting up %.0fM-row out-of-core benchmark (DuckDB %d threads, %s cap) ===\n",
            N_ROWS / 1e6, THREADS, MEM_LIM))

con_duck <- dbConnect(duckdb::duckdb(), dbdir = duck_file)
on.exit(dbDisconnect(con_duck, shutdown = TRUE), add = TRUE)
dbExecute(con_duck, sprintf("PRAGMA threads=%d", THREADS))
dbExecute(con_duck, sprintf("SET memory_limit='%s'", MEM_LIM))
dbExecute(con_duck, sprintf("SET temp_directory='%s'", duck_tmp))

t0 <- proc.time()
dbExecute(con_duck, paste("CREATE TABLE bench_100m AS", table_100m_sql(N_ROWS)))
cat(sprintf("DuckDB materialized in %.2fs (on-disk .duckdb: %.2f GiB)\n",
            (proc.time() - t0)[3], file.info(duck_file)$size / (1024^3)))

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
  dbExecute(con_dbx, paste("CREATE TEMP TABLE bench_100m AS", table_100m_sql(N_ROWS)))
  cat(sprintf("Databricks materialized in %.2fs\n\n", (proc.time() - t0)[3]))
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
  list(name = "01_OLS_HC1_100M", has_upstream = TRUE,
       fn = function(con) dbreg(yg ~ x1 + x2 + x3, conn = con, table = "bench_100m", vcov = "hc1")),
  list(name = "02_1FE_10k_Cluster_100M", has_upstream = TRUE,
       fn = function(con) dbreg(yg ~ x1 + x2 + x3 | teacher, conn = con, table = "bench_100m", vcov = ~school)),
  list(name = "03_Compress_OLS_2FE_100M", has_upstream = TRUE,
       fn = function(con) dbreg(permanent_wastage ~ FTPTWK + SEX + AGE | district + YEAR_D, conn = con, table = "bench_100m")),
  list(name = "04_Binscatter_NoFE_100M", has_upstream = TRUE,
       fn = function(con) dbbinsreg(yg ~ x1 + x2, conn = con, table = "bench_100m", nbins = 20, plot = FALSE)),
  list(name = "05_Binscatter_1FE_100M", has_upstream = TRUE,
       fn = function(con) dbbinsreg(yg ~ x1 + year | school, conn = con, table = "bench_100m", nbins = 20, plot = FALSE)),
  list(name = "06_2FE_10kx100_100M", has_upstream = TRUE,
       fn = function(con) dbreg(yg ~ x1 + x2 + x3 | teacher + school, conn = con, table = "bench_100m", vcov = "hc1")),
  list(name = "07_GLM_Logit_Compress_100M", has_upstream = FALSE,
       fn = function(con) dbglm(permanent_wastage ~ FTPTWK + SEX + AGE | district + YEAR_D, family = binomial(), conn = con, table = "bench_100m"))
)

run_single <- function(ver, engine, con, fn) {
  use_version(ver)
  t0 <- proc.time()
  fit <- tryCatch(fn(con), error = function(e) structure(list(error = conditionMessage(e)), class = "bench_err"))
  elapsed <- (proc.time() - t0)[3]
  if (inherits(fit, "bench_err")) {
    cat(sprintf("  [%-8s | %-10s] ERROR (%.2fs): %s\n", ver, engine, elapsed, fit$error))
    return(list(elapsed = NA_real_, fit = NULL, nobs = NA_real_))
  }
  n_obs <- if (inherits(fit, "dbbinsreg")) as.numeric(fit$model$nobs_orig) else as.numeric(nobs(fit))
  cat(sprintf("  [%-8s | %-10s] %.3fs (nobs=%.0f)\n", ver, engine, elapsed, n_obs))
  list(elapsed = elapsed, fit = fit, nobs = n_obs)
}

extract_coefs <- function(fit) {
  if (is.null(fit)) return(NULL)
  if (inherits(fit, "dbbinsreg")) return(as.numeric(fit$points$fit))
  as.numeric(coef(fit))
}

all_rows <- list()
for (c_info in cases) {
  cat(sprintf("--- Running 12GB Case: %s ---\n", c_info$name))
  f_duck <- run_single("fork", "DuckDB", con_duck, c_info$fn)
  u_duck <- if (has_upstream && c_info$has_upstream) run_single("upstream", "DuckDB", con_duck, c_info$fn) else list(elapsed = NA_real_, fit = NULL, nobs = NA_real_)
  f_dbx  <- if (!is.null(con_dbx)) run_single("fork", "Databricks", con_dbx, c_info$fn) else list(elapsed = NA_real_, fit = NULL, nobs = NA_real_)
  u_dbx  <- if (!is.null(con_dbx) && has_upstream && c_info$has_upstream) run_single("upstream", "Databricks", con_dbx, c_info$fn) else list(elapsed = NA_real_, fit = NULL, nobs = NA_real_)
  c_fd <- extract_coefs(f_duck$fit); c_fb <- extract_coefs(f_dbx$fit)
  coef_diff <- if (!is.null(c_fd) && !is.null(c_fb) && length(c_fd) == length(c_fb)) max(abs(c_fd - c_fb)) else NA_real_

  all_rows[[length(all_rows) + 1L]] <- data.frame(
    case = c_info$name,
    fork_duck_s = round(f_duck$elapsed, 3),
    upst_duck_s = round(u_duck$elapsed, 3),
    fork_dbx_s  = round(f_dbx$elapsed, 3),
    upst_dbx_s  = round(u_dbx$elapsed, 3),
    fork_nobs_duck = f_duck$nobs,
    upst_nobs_duck = u_duck$nobs,
    max_coef_diff_duck_dbx = coef_diff,
    stringsAsFactors = FALSE
  )
  write.csv(do.call(rbind, all_rows), CSV_OUT, row.names = FALSE)
}

cat("\n=== DONE ===\n")
print(do.call(rbind, all_rows), row.names = FALSE)
