#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(fixest)
  library(odbc)
  library(pkgload)
})

script_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
repo = if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."))
} else {
  normalizePath(".")
}

warehouse_id = Sys.getenv("DATABRICKS_WAREHOUSE_ID")
http_path = if (nzchar(warehouse_id)) {
  paste0("/sql/1.0/warehouses/", warehouse_id)
} else {
  Sys.getenv("DATABRICKS_SQL_PATH")
}
if (!nzchar(http_path)) {
  stop("Set DATABRICKS_WAREHOUSE_ID or DATABRICKS_SQL_PATH.")
}

schema = Sys.getenv("DBREG_DATABRICKS_BENCHMARK_SCHEMA")
if (!nzchar(schema)) {
  stop("Set DBREG_DATABRICKS_BENCHMARK_SCHEMA to a writable catalog.schema.")
}

n_rows = as.numeric(Sys.getenv("DBREG_BENCHMARK_ROWS", "100000"))
if (!is.finite(n_rows) || n_rows < 1000 || n_rows != floor(n_rows)) {
  stop("DBREG_BENCHMARK_ROWS must be an integer of at least 1000.")
}

pkgload::load_all(repo, quiet = TRUE)
con = dbConnect(
  odbc::databricks(),
  httpPath = http_path,
  useNativeQuery = TRUE,
  DefaultStringColumnLength = 65536
)
on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)

table = paste0(schema, ".dbreg_bounded_ap_benchmark")
dbExecute(con, sprintf(
  paste0(
    "CREATE OR REPLACE TABLE %s USING DELTA AS SELECT ",
    "CAST(id AS BIGINT) AS row_id, ",
    "CAST(id %% 100 AS BIGINT) AS fe_1, ",
    "CAST(FLOOR(id / 13) %% 40 AS BIGINT) AS fe_2, ",
    "CAST(FLOOR(id / 37) %% 20 AS BIGINT) AS fe_3, ",
    "RAND(51) AS x_1, RAND(52) AS x_2, ",
    "CAST(0.7 * RAND(51) - 0.4 * RAND(52) + ",
    "SIN((id %% 100) * 0.1) + COS((FLOOR(id / 100) %% 20) * 0.2) + ",
    "SIN((FLOOR(id / 17) %% 113) * 0.03) + ",
    "RANDN(53) AS DOUBLE) AS y FROM RANGE(0, %.0f)"
  ),
  table,
  n_rows
))

formula = y ~ x_1 + x_2 | fe_1 + fe_2 + fe_3
reference_data = dbGetQuery(
  con,
  paste("SELECT y, x_1, x_2, fe_1, fe_2, fe_3 FROM", table)
)
reference = fixest::feols(formula, data = reference_data, vcov = "iid", fixef.rm = "none")

options(
  dbreg.ap_tol = 1e-7,
  dbreg.ap_max_iter = 100L,
  dbreg.spark_ap_chunk_cycles = 3L
)
started = proc.time()[["elapsed"]]
fit = dbreg::dbreg(
  formula,
  conn = con,
  table = table,
  strategy = "demean",
  vcov = "iid",
  fixef_rm = "none"
)
elapsed = proc.time()[["elapsed"]] - started
coefficient_error = max(abs(coef(fit) - coef(reference)))

result = data.frame(
  git_sha = system2("git", c("-C", repo, "rev-parse", "--short", "HEAD"), stdout = TRUE),
  rows = n_rows,
  elapsed_seconds = elapsed,
  ap_cycles = fit$ap_iterations,
  ap_max_abs = fit$ap_max_abs,
  coefficient_error = coefficient_error,
  stringsAsFactors = FALSE
)
print(result, row.names = FALSE)
cat("Fixture table:", table, "\n")

if (!is.finite(coefficient_error) || coefficient_error > 1e-6) {
  stop("Coefficient agreement with fixest exceeded 1e-6.", call. = FALSE)
}
