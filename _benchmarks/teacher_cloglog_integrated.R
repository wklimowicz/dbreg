# Benchmark the actual package implementation: no namespace substitutions and
# the user's original uncentered formula. Reuse the exact strict oracle/sample.
# Rscript _benchmarks/teacher_cloglog_integrated.R /tmp/teacher-cloglog-5m /tmp/teacher-cloglog-integrated 2
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("_benchmarks/teacher_cloglog_profile.R")
args = commandArgs(TRUE)
input = args[[1]]; out = args[[2]]
repetitions = if (length(args) > 2L) as.integer(args[[3]]) else 2L
dir.create(file.path(out, "sql"), recursive = TRUE, showWarnings = FALSE)
reference = readRDS(file.path(input, "fixest_strict_1.rds"))
summary = read.csv(file.path(input, "data_summary.csv"))
log = new.env(parent = emptyenv())
log[["out"]] = out; log[["id"]] = 0L; log[["tag"]] = paste0("teacher_integrated_", Sys.getpid())
log[["run"]] = "setup"; log[["queries"]] = list()
run = function() {
  native = DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
  on.exit(DBI::dbDisconnect(native))
  remote = methods::new("teacher_profile_connection", native = native, log = log)
  DBI::dbExecute(native, "SET use_cached_result=false")
  local = DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB", temp_directory = file.path(out, "duckdb_tmp"))))
  on.exit(DBI::dbDisconnect(local, shutdown = TRUE), add = TRUE)
  settings = DBI::dbGetQuery(local, "SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")
  stopifnot(as.integer(settings[["threads"]]) == 4L)
  write.csv(settings, file.path(out, "duckdb_settings.csv"), row.names = FALSE)
  table = paste0("teacher_integrated_input_", Sys.getpid())
  sql = teacher_cloglog_sql(summary[["n"]])
  stopifnot(identical(sql, paste(readLines(file.path(input, "data_generation.sql")), collapse = "\n")))
  DBI::dbExecute(native, paste("CREATE TEMP TABLE", table, "AS", sql))
  counts = DBI::dbGetQuery(native, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits FROM", table))
  stopifnot(as.numeric(counts[["n"]]) == summary[["n"]], as.numeric(counts[["exits"]]) == summary[["exits"]])
  DBI::dbExecute(local, paste0("CREATE TEMP TABLE teachers AS SELECT * FROM read_parquet(",
    DBI::dbQuoteString(local, file.path(input, "teachers.parquet")), ")"))
  formula = wastage ~ log_outside_salary | spine_band_year + district + Experience + AGE + SEX + FTPTWK
  environment(formula) = baseenv()
  results = list()
  for (engine in c("Databricks", "DuckDB")) for (rep in seq_len(repetitions)) {
    log[["run"]] = paste(engine, "integrated", rep, sep = "_")
    conn = if (engine == "Databricks") remote else local
    source = if (engine == "Databricks") table else "teachers"
    DBI::dbGetQuery(if (engine == "Databricks") native else local, "SELECT 1 AS warmup")
    gc()
    cat("START", log[["run"]], format(Sys.time()), "\n"); flush.console()
    timing = system.time({fit = dbglm(formula, conn = conn, table = source,
      family = binomial("cloglog"), vcov = ~district, verbose = TRUE)})
    row = data.frame(engine = engine, repetition = rep, elapsed = timing[["elapsed"]],
      coefficient = unname(coef(fit)[[1]]), se = sqrt(vcov(fit)[1, 1]),
      coef_error = max(abs(coef(fit) - reference[["coefficients"]])),
      se_error = max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(reference[["vcov"]])))),
      nobs = nobs(fit), iterations = fit[["iterations"]], score = fit[["score_max"]],
      halvings = fit[["step_halvings"]], method = fit[["projection_method"]],
      deviance = deviance(fit), design_center = unname(fit[["design_center"]][[1]]))
    results[[length(results) + 1L]] = row
    write.csv(do.call(rbind, results), file.path(out, "fits.csv"), row.names = FALSE)
    saveRDS(fit, file.path(out, paste0(log[["run"]], ".rds")))
    print(row); flush.console()
    stopifnot(nobs(fit) == summary[["n"]], row[["coef_error"]] < 1e-8,
      row[["se_error"]] < 1e-7, abs(deviance(fit) - reference[["deviance"]]) < 1e-5,
      fit[["score_max"]] <= fit[["control"]][["score_tol"]])
  }
  DBI::dbExecute(native, paste("DROP TEMP TABLE", table))
}
run()
