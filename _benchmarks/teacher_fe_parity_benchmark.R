# Ordinary FE regression on the exact 5M teacher-shaped GLM benchmark input.
# Rscript _benchmarks/teacher_fe_parity_benchmark.R INPUT OUT local|remote [REPS] [MAX_ITER]
# Run local and remote sequentially, without competing tests/benchmarks.
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
args = commandArgs(TRUE)
input = args[[1]]
out = args[[2]]
mode = match.arg(args[[3]], c("local", "remote"))
repetitions = if (length(args) > 3L) as.integer(args[[4]]) else 2L
max_iter = if (length(args) > 4L) as.integer(args[[5]]) else 1000L
dir.create(out, recursive = TRUE, showWarnings = FALSE)

run = function() {
  saved = options(dbreg.ap_tol = 1e-7, dbreg.ap_max_iter = max_iter)
  on.exit(options(saved))
  old_threads = fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(4)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)
  formula = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(formula) = baseenv()
  summary = read.csv(file.path(input, "data_summary.csv"))
  results = list()
  reference_path = file.path(out, "feols_reference.rds")
  if (mode == "local") {
    conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB",
      temp_directory = file.path(out, "duckdb_tmp"))))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(conn, paste0("CREATE TEMP TABLE teachers AS SELECT * FROM read_parquet(",
      DBI::dbQuoteString(conn, file.path(input, "teachers.parquet")), ")"))
    d = DBI::dbGetQuery(conn, "SELECT * FROM teachers")
    stopifnot(nrow(d) == summary[["n"]], sum(d[["wastage"]]) == summary[["exits"]])
    reference = fixest::feols(formula, d, vcov = ~district, fixef.rm = "none",
      fixef.tol = 1e-10, fixef.iter = 10000, ssc = fixest::ssc(K.fixef = "full"), nthreads = 4)
    saveRDS(list(coef = coef(reference), vcov = vcov(reference), nobs = nobs(reference)), reference_path)
    engines = c("fixest", "DuckDB")
    table = "teachers"
    metadata = c(capture.output(sessionInfo()), "DuckDB threads: 4; fixest threads: 4",
      paste("dbreg AP tolerance: 1e-7; max iterations:", max_iter),
      "fixest timed fits: package-default demeaning tolerance; reference: 1e-10")
    writeLines(metadata, file.path(out, "local_environment.txt"))
  } else {
    conn = DBI::dbConnect(odbc::databricks(),
      driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
      httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
    on.exit(DBI::dbDisconnect(conn), add = TRUE)
    DBI::dbExecute(conn, "SET use_cached_result=false")
    table = paste0("teacher_fe_parity_", Sys.getpid())
    sql = teacher_cloglog_sql(summary[["n"]])
    stopifnot(identical(sql, paste(readLines(file.path(input, "data_generation.sql")), collapse = "\n")))
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", sql))
    counts = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits FROM", table))
    stopifnot(as.numeric(counts[["n"]]) == summary[["n"]], as.numeric(counts[["exits"]]) == summary[["exits"]])
    engines = "Databricks"
    writeLines(c(capture.output(sessionInfo()),
      capture.output(DBI::dbGetQuery(conn, "SELECT current_version() AS version")),
      "Result reuse disabled; generation and warehouse startup excluded."), file.path(out, "remote_environment.txt"))
  }
  ref = readRDS(reference_path)
  for (rep in seq_len(repetitions)) for (engine in if (rep %% 2L) engines else rev(engines)) {
    DBI::dbGetQuery(conn, "SELECT 1 AS warmup")
    gc()
    message("START ", engine, " ", rep, " ", Sys.time())
    elapsed = system.time({attempt = tryCatch({
      fit = if (engine == "fixest") {
        fixest::feols(formula, d, vcov = ~district, fixef.rm = "none",
          ssc = fixest::ssc(K.fixef = "full"), nthreads = 4)
      } else dbreg(formula, conn = conn, table = table, strategy = "demean", vcov = ~district, verbose = TRUE)
      list(fit = fit, error = "")
    }, error = function(e) list(fit = NULL, error = conditionMessage(e)))})[["elapsed"]]
    fit = attempt[["fit"]]
    if (is.null(fit)) {
      # A capped diagnostic is an unsuccessful fit, never a valid fit timing.
      row = data.frame(engine = engine, repetition = rep, elapsed = elapsed, nobs = NA_real_,
        coefficient = NA_real_, se = NA_real_, coef_error = NA_real_, se_error = NA_real_,
        vcov_error = NA_real_, projection_iterations = NA_integer_,
        status = if (grepl("did not converge", attempt[["error"]], fixed = TRUE)) "not_converged" else "error",
        max_iter = max_iter, error = attempt[["error"]])
      results[[length(results) + 1L]] = row
      write.csv(do.call(rbind, results), file.path(out, paste0(mode, "_fits.csv")), row.names = FALSE)
      print(row)
      next
    }
    beta = coef(fit)
    V = vcov(fit)
    row = data.frame(engine = engine, repetition = rep, elapsed = elapsed, nobs = nobs(fit),
      coefficient = unname(beta[[1]]), se = sqrt(V[1, 1]),
      coef_error = max(abs(beta - ref[["coef"]])),
      se_error = max(abs(sqrt(diag(V)) - sqrt(diag(ref[["vcov"]])))),
      vcov_error = max(abs(V - ref[["vcov"]])),
      projection_iterations = if (is.null(fit[["ap_iterations"]])) NA_integer_ else fit[["ap_iterations"]],
      status = "converged", max_iter = max_iter, error = "")
    results[[length(results) + 1L]] = row
    write.csv(do.call(rbind, results), file.path(out, paste0(mode, "_fits.csv")), row.names = FALSE)
    print(row)
    stopifnot(identical(names(beta), names(ref[["coef"]])), nobs(fit) == ref[["nobs"]],
      all(abs(beta - ref[["coef"]]) <= 1e-7 + 1e-6 * abs(ref[["coef"]])),
      all(abs(sqrt(diag(V)) - sqrt(diag(ref[["vcov"]]))) <= 1e-8 + 1e-5 * sqrt(diag(ref[["vcov"]]))),
      all(abs(V - ref[["vcov"]]) <= 1e-8 + 1e-5 * abs(ref[["vcov"]])))
  }
  if (mode == "remote") DBI::dbExecute(conn, paste("DROP TEMP TABLE", table))
}
run()
