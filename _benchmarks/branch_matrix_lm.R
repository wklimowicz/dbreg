# Synthetic 5M-row, ten-slope LM counterpart to teacher_cloglog_ten_slopes.R.
# Run engines sequentially. Never edit a running copy of this script.
# Rscript _benchmarks/branch_matrix_lm.R REPO INPUT OUT ENGINE LABEL SHA REPS MAX_ITER [MAX_CELLS]
# ENGINE: local (strict reference + fixest + DuckDB), fixest, duckdb, databricks.
args = commandArgs(TRUE)
stopifnot(length(args) %in% c(8L, 9L))
repo = normalizePath(args[[1]])
input = normalizePath(args[[2]])
out = args[[3]]
engine = match.arg(args[[4]], c("local", "fixest", "duckdb", "databricks"))
label = args[[5]]
sha = args[[6]]
repetitions = as.integer(args[[7]])
max_iter = as.integer(args[[8]])
max_cells = if (length(args) == 9L) as.numeric(args[[9]]) else 8e6
stopifnot(repetitions > 0L, max_iter > 0L)
pkgload::load_all(repo)
source("_benchmarks/teacher_cloglog_data.R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

run = function() {
  saved = options(dbreg.ap_tol = 1e-7, dbreg.ap_max_iter = max_iter,
    dbreg.ap_moment_max_cells = max_cells)
  on.exit(options(saved))
  fixest::setFixest_nthreads(4)
  n = 5000000L
  regressors = c(paste0("outside_", 1:5), paste0("teacher_", 1:5))
  formula = as.formula(paste("wastage ~", paste(regressors, collapse = " + "),
    "| spine_band_year + district + Experience + AGE + SEX + FTPTWK"))
  environment(formula) = baseenv()
  remote = engine == "databricks"
  spill = file.path(out, paste0(label, "_spill"))
  if (!remote) dir.create(spill, recursive = TRUE, showWarnings = FALSE)
  conn = if (remote) DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE,
    DefaultStringColumnLength = 65536) else
    DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB", temp_directory = spill)))
  on.exit(if (remote) DBI::dbDisconnect(conn) else DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  if (!remote) on.exit(unlink(spill, recursive = TRUE), add = TRUE)
  if (remote) DBI::dbExecute(conn, "SET use_cached_result=false")
  table = paste0("branch_matrix_", Sys.getpid())
  extra = paste0(teacher_cloglog_salary_sql(), ",",
    "CASE WHEN Experience<5 THEN 1 WHEN Experience<10 THEN 2 WHEN Experience<20 THEN 3 ",
    "WHEN Experience<30 THEN 4 ELSE 5 END AS experience_group")
  original = if (remote) paste0("(", teacher_cloglog_sql(n), ") original") else
    paste0("read_parquet(", DBI::dbQuoteString(conn, file.path(input, "teachers.parquet")), ")")
  # Explicit numeric columns give every branch exactly the same design matrix;
  # dbreg's formula validator does not accept factor() calls in interactions.
  wage_terms = c(sprintf("CASE WHEN experience_group=%d THEN log_outside_salary ELSE 0.0 END AS outside_%d", 1:5, 1:5),
    sprintf("CASE WHEN experience_group=%d THEN log_teacher_salary ELSE 0.0 END AS teacher_%d", 1:5, 1:5))
  DBI::dbExecute(conn, paste0("CREATE TEMP TABLE ", table, " AS SELECT *,", paste(wage_terms, collapse = ","),
    " FROM (SELECT *,", extra, " FROM ", original, ") expanded"))
  drop = if (remote) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
  on.exit(DBI::dbExecute(conn, paste(drop, table)), add = TRUE, after = FALSE)
  summary = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits,",
    "COUNT(DISTINCT id) AS ids,MIN(id) AS lo,MAX(id) AS hi,",
    "COUNT(DISTINCT experience_group) AS groups,COUNT(DISTINCT Experience) AS experience_levels,",
    "AVG(log_outside_salary) AS outside_mean,AVG(log_teacher_salary) AS teacher_mean FROM", table))
  for (name in names(summary)) summary[[name]] = as.numeric(summary[[name]])
  stopifnot(summary[["n"]] == n, summary[["exits"]] == 105021,
    summary[["ids"]] == n, summary[["lo"]] == 0, summary[["hi"]] == n - 1L,
    summary[["groups"]] == 5, summary[["experience_levels"]] == 44)
  snapshot = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "WHERE id % 49999=0 ORDER BY id"))
  for (name in names(snapshot)) if (inherits(snapshot[[name]], "integer64")) snapshot[[name]] = as.numeric(snapshot[[name]])
  if (engine == "local") {
    saveRDS(snapshot, file.path(out, "input_snapshot.rds"))
    write.csv(summary, file.path(out, "input_summary.csv"), row.names = FALSE)
  } else {
    ref_summary = read.csv(file.path(out, "input_summary.csv"))
    for (name in names(summary)) stopifnot(abs(summary[[name]] - ref_summary[[name]]) <=
      if (grepl("_mean$", name)) 1e-10 else 0)
    stopifnot(isTRUE(all.equal(snapshot, readRDS(file.path(out, "input_snapshot.rds")),
      check.attributes = FALSE, tolerance = 1e-14)))
  }
  write.csv(summary, file.path(out, paste0(label, "_input.csv")), row.names = FALSE)
  metadata = c(paste("Source:", sha, "label:", label, "repo:", repo),
    paste("AP tolerance: 1e-7; max iterations:", max_iter),
    paste("Bounded FE rectangular transfer-cell limit:", max_cells),
    "Synthetic only. 5M rows, ten slopes, six FEs, 296 district clusters; full FE parameter counting.",
    "Timers include input processing, fitting and covariance; exclude data preparation, connections and reference fits.",
    "Four local threads; BLAS=1. Native result reuse disabled; other caches may be warm.",
    capture.output(sessionInfo()),
    if (remote) capture.output(DBI::dbGetQuery(conn, "SELECT current_version() AS version")) else
      capture.output(DBI::dbGetQuery(conn, "SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")))
  writeLines(metadata, file.path(out, paste0(label, "_environment.txt")))
  if (engine %in% c("local", "fixest")) {
    data = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "ORDER BY id"))
    for (group in 1:5) {
      stopifnot(identical(data[[paste0("outside_", group)]],
        data[["log_outside_salary"]] * (data[["experience_group"]] == group)),
        identical(data[[paste0("teacher_", group)]],
        data[["log_teacher_salary"]] * (data[["experience_group"]] == group)))
    }
    if (engine == "local") {
      message("Building strict LM reference")
      strict = fixest::feols(formula, data, vcov = ~district, fixef.rm = "none",
        ssc = fixest::ssc(K.fixef = "full"), fixef.tol = 1e-10, fixef.iter = 10000)
      stopifnot(length(coef(strict)) == 10L, nobs(strict) == n)
      saveRDS(list(coef = coef(strict), vcov = vcov(strict), nobs = nobs(strict)), file.path(out, "reference.rds"))
    }
  }
  ref = readRDS(file.path(out, "reference.rds"))
  engines = if (engine == "local") c("fixest", "duckdb") else engine
  rows = list()
  for (rep in seq_len(repetitions)) for (current in if (rep %% 2L) engines else rev(engines)) {
    DBI::dbGetQuery(conn, "SELECT 1 AS warmup")
    gc()
    message("START ", label, " ", current, " ", rep, " ", Sys.time())
    elapsed = system.time({attempt = tryCatch({
      fit = if (current == "fixest") fixest::feols(formula, data, vcov = ~district,
        fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full")) else
        dbreg(formula, conn = conn, table = table, strategy = "demean", vcov = ~district, verbose = TRUE)
      list(fit = fit, error = "")
    }, error = function(e) list(fit = NULL, error = conditionMessage(e)))})[["elapsed"]]
    fit = attempt[["fit"]]
    row = data.frame(label = label, sha = if (current == "fixest") as.character(packageVersion("fixest")) else sha,
      engine = current, repetition = rep, elapsed = elapsed,
      max_iter = max_iter, max_cells = max_cells, nobs = NA_real_, slopes = NA_integer_, projection_iterations = NA_integer_,
      coef_error = NA_real_, vcov_error = NA_real_, se_error = NA_real_, parity = FALSE,
      status = if (grepl("did not converge", attempt[["error"]], fixed = TRUE)) "not_converged" else "error",
      error = attempt[["error"]])
    if (!is.null(fit)) {
      beta = coef(fit)
      V = vcov(fit)
      same_names = identical(names(beta), names(ref[["coef"]]))
      row[["nobs"]] = as.numeric(nobs(fit))
      row[["slopes"]] = length(beta)
      row[["projection_iterations"]] = if (is.null(fit[["ap_iterations"]])) NA_integer_ else fit[["ap_iterations"]]
      if (same_names) {
        row[["coef_error"]] = max(abs(beta - ref[["coef"]]))
        row[["vcov_error"]] = max(abs(V - ref[["vcov"]]))
        row[["se_error"]] = max(abs(sqrt(diag(V)) - sqrt(diag(ref[["vcov"]]))))
      }
      row[["parity"]] = same_names && row[["nobs"]] == n && row[["slopes"]] == 10L &&
        all(abs(beta - ref[["coef"]]) <= 1e-7 + 1e-6 * abs(ref[["coef"]])) &&
        all(abs(V - ref[["vcov"]]) <= 1e-8 + 1e-5 * abs(ref[["vcov"]])) &&
        all(abs(sqrt(diag(V)) - sqrt(diag(ref[["vcov"]]))) <= 1e-8 + 1e-5 * sqrt(diag(ref[["vcov"]])))
      row[["status"]] = if (isTRUE(row[["parity"]])) "converged" else "parity_failed"
      saveRDS(list(coef = beta, vcov = V, nobs = row[["nobs"]]),
        file.path(out, paste0(label, "_", current, "_", rep, ".rds")))
    }
    rows[[length(rows) + 1L]] = row
    write.csv(do.call(rbind, rows), file.path(out, paste0(label, "_fits.csv")), row.names = FALSE)
    print(row)
  }
}
run()
