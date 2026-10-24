# Sequential, synthetic-only 5M-row / ten-slope cloglog comparison.
# Rscript _benchmarks/teacher_cloglog_ten_slopes.R INPUT OUT MODE [REPS] [N]
# MODE: duckdb, fixest, databricks, reference, sparse, compare.
# Benchmark engines use their ordinary defaults. Strict references are separate.
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("inst/tinytest/helper_glm_reference.R")
source("_benchmarks/teacher_cloglog_sparse_reference.R")
args = commandArgs(TRUE)
input = args[[1]]
out = args[[2]]
mode = match.arg(args[[3]], c("duckdb", "fixest", "databricks", "reference", "sparse", "compare"))
repetitions = if (length(args) > 3L) as.integer(args[[4]]) else 2L
n = if (length(args) > 4L) as.integer(args[[5]]) else 5000000L
stopifnot(repetitions > 0L, n >= 10000L, n <= 5000000L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
formula = wastage ~ log_outside_salary:factor(experience_group) +
  log_teacher_salary:factor(experience_group) |
  spine_band_year + district + Experience + AGE + SEX + FTPTWK
environment(formula) = baseenv()
fes = c("spine_band_year", "district", "Experience", "AGE", "SEX", "FTPTWK")
extra = paste0(teacher_cloglog_salary_sql(), ",",
  "CASE WHEN Experience<5 THEN 1 WHEN Experience<10 THEN 2 WHEN Experience<20 THEN 3 ",
  "WHEN Experience<30 THEN 4 ELSE 5 END AS experience_group")

summarize_fixest = function(fit, data, elapsed, repetition = NA_integer_) {
  list(engine = "fixest", repetition = repetition, elapsed = elapsed,
    coef = coef(fit), vcov = vcov(fit), logLik = as.numeric(logLik(fit)), nobs = nobs(fit),
    n_input = nrow(data), exits = sum(data[["wastage"]]),
    keys = as.numeric(data[["id"]][fixest::obs(fit)]), model_df = fit[["nparams"]],
    covariance_df = attr(vcov(fit, attr = TRUE), "df.K"), iterations = fit[["iterations"]],
    converged = isTRUE(fit[["convStatus"]]), score = NA_real_, projection = NA_real_,
    refinements = attr(fit, "reference_refinements"))
}

compare_saved = function() {
  references = c(fixest = "reference.rds", sparse = "sparse_reference.rds")
  results = list()
  checks = list()
  for (label in names(references)) {
    path = file.path(out, references[[label]])
    if (!file.exists(path)) next
    ref = readRDS(path)
    files = list.files(out, "^(duckdb|fixest|databricks)_[0-9]+[.]rds$", full.names = TRUE)
    strict_file = file.path(out, "reference.rds")
    if (label == "sparse" && file.exists(strict_file)) files = c(files, strict_file)
    for (file in files) {
      fit = readRDS(file)
      if (basename(file) == "reference.rds") {
        fit[["engine"]] = "fixest_strict"
        fit[["repetition"]] = 0L
      }
      stopifnot(fit[["n_input"]] == ref[["n_input"]], fit[["exits"]] == ref[["exits"]])
      test = c(names = identical(names(fit[["coef"]]), names(ref[["coef"]])),
        nobs = fit[["nobs"]] == ref[["nobs"]], sample = identical(fit[["keys"]], sort(ref[["keys"]])),
        coef = all(abs(fit[["coef"]] - ref[["coef"]]) <= 1e-7 + 1e-6 * abs(ref[["coef"]])),
        vcov = all(abs(fit[["vcov"]] - ref[["vcov"]]) <= 1e-8 + 1e-5 * abs(ref[["vcov"]])),
        se = all(abs(sqrt(diag(fit[["vcov"]])) - sqrt(diag(ref[["vcov"]]))) <=
          1e-8 + 1e-6 * sqrt(diag(ref[["vcov"]]))),
        likelihood = abs(fit[["logLik"]] - ref[["logLik"]]) / fit[["nobs"]] < 1e-9,
        model_df = fit[["model_df"]] == ref[["model_df"]],
        covariance_df = fit[["covariance_df"]] == ref[["covariance_df"]],
        converged = fit[["converged"]])
      if (fit[["engine"]] %in% c("duckdb", "databricks")) test = c(test,
        score = fit[["score"]] <= fit[["score_tol"]],
        projection = fit[["projection"]] <= fit[["fe_tol"]])
      checks[[length(checks) + 1L]] = data.frame(engine = fit[["engine"]],
        repetition = fit[["repetition"]], reference = label, check = names(test), passed = unname(test))
      results[[length(results) + 1L]] = data.frame(engine = fit[["engine"]],
        repetition = fit[["repetition"]], reference = label, nobs = fit[["nobs"]],
        coefficients = length(fit[["coef"]]), elapsed = fit[["elapsed"]],
        iterations = fit[["iterations"]], coef_error = max(abs(fit[["coef"]] - ref[["coef"]])),
        vcov_error = max(abs(fit[["vcov"]] - ref[["vcov"]])),
        se_error = max(abs(sqrt(diag(fit[["vcov"]])) - sqrt(diag(ref[["vcov"]])))),
        likelihood_error = abs(fit[["logLik"]] - ref[["logLik"]]), passed = all(test))
    }
  }
  stopifnot(length(results) > 0L)
  write.csv(do.call(rbind, checks), file.path(out, "checks.csv"), row.names = FALSE)
  result = do.call(rbind, results)
  write.csv(result, file.path(out, "comparison.csv"), row.names = FALSE)
  print(result)
  if (any(!result[["passed"]])) stop("Some parity gates failed; inspect checks.csv and comparison.csv.")
  invisible(result)
}

run = function() {
  if (mode == "compare") return(compare_saved())
  fixest::setFixest_nthreads(4)
  remote = mode == "databricks"
  spill = file.path(out, "spill", paste0(mode, "-", Sys.getpid()))
  if (!remote) dir.create(spill, recursive = TRUE, showWarnings = FALSE)
  conn = if (remote) DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE,
    DefaultStringColumnLength = 65536) else
    DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB", temp_directory = spill)))
  on.exit(if (remote) DBI::dbDisconnect(conn) else DBI::dbDisconnect(conn, shutdown = TRUE))
  if (!remote) on.exit(unlink(spill, recursive = TRUE), add = TRUE)
  if (remote) DBI::dbExecute(conn, "SET use_cached_result=false")
  table = paste0("teacher_ten_", Sys.getpid())
  retained = paste0(table, "_retained")
  original = if (remote) paste0("(", teacher_cloglog_sql(n), ") original") else
    paste0("read_parquet(", DBI::dbQuoteString(conn, file.path(input, "teachers.parquet")), ") WHERE id<", n)
  DBI::dbExecute(conn, paste0("CREATE TEMP TABLE ", table, " AS SELECT *,", extra, " FROM ", original))
  drop = if (remote) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
  on.exit(for (name in c(retained, table)) DBI::dbExecute(conn, paste(drop, name)), add = TRUE, after = FALSE)
  summary = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits,",
    "COUNT(DISTINCT id) AS ids,MIN(id) AS lo,MAX(id) AS hi,",
    "COUNT(DISTINCT experience_group) AS groups,COUNT(DISTINCT Experience) AS experience_levels,",
    "AVG(log_outside_salary) AS outside_mean,AVG(log_teacher_salary) AS teacher_mean FROM", table))
  stopifnot(as.numeric(summary[["n"]]) == n, as.numeric(summary[["ids"]]) == n,
    as.numeric(summary[["lo"]]) == 0, as.numeric(summary[["hi"]]) == n - 1L,
    as.numeric(summary[["groups"]]) == 5L, as.numeric(summary[["experience_levels"]]) == 44L)
  if (n == 5000000L) stopifnot(as.numeric(summary[["exits"]]) == 105021L)
  write.csv(summary, file.path(out, paste0(mode, "_input.csv")), row.names = FALSE)
  snapshot = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "WHERE id % 49999=0 ORDER BY id"))
  for (name in names(snapshot)) if (inherits(snapshot[[name]], "integer64")) snapshot[[name]] = as.numeric(snapshot[[name]])
  if (remote) {
    # Fail before starting any timer if the backend expressions generated
    # different inputs. Explicit DOUBLE casts above avoid Spark DECIMAL division.
    local_summary = read.csv(file.path(out, "duckdb_input.csv"))
    for (name in names(summary)) {
      tolerance = if (grepl("_mean$", name)) 1e-10 else 0
      stopifnot(abs(as.numeric(summary[[name]]) - local_summary[[name]]) <= tolerance)
    }
    local_snapshot = readRDS(file.path(out, "duckdb_input_snapshot.rds"))
    stopifnot(isTRUE(all.equal(snapshot, local_snapshot, check.attributes = FALSE, tolerance = 1e-14)))
    message("Native input matches local counts, wage means and all sampled row values")
  }
  saveRDS(snapshot, file.path(out, paste0(mode, "_input_snapshot.rds")))
  environment = c(capture.output(sessionInfo()), paste("mode:", mode, "rows:", n),
    "Synthetic only. Two wages x five experience groups (0-4, 5-9, 10-19, 20-29, 30-43).",
    "Six original FEs, including 44 single-year experience levels; district-clustered covariance.",
    "DuckDB and fixest threads=4; BLAS=1. Run engines sequentially with data already prepared.",
    "Timers include input processing, fitting, covariance and dbglm sample export; exclude source loading and reference fits.",
    "Default estimator tolerances/initialization in benchmark modes. No namespace overrides.",
    if (remote) capture.output(DBI::dbGetQuery(conn, "SELECT current_version() AS version")) else
      capture.output(DBI::dbGetQuery(conn, "SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")))
  writeLines(environment, file.path(out, paste0(mode, "_environment.txt")))
  if (mode %in% c("fixest", "reference", "sparse")) {
    data = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "ORDER BY id"))
    stopifnot(identical(as.numeric(data[["id"]]), as.numeric(seq_len(n) - 1L)))
    if (mode %in% c("reference", "sparse")) {
      for (name in c("log_outside_salary", "log_teacher_salary")) data[[name]] = data[[name]] - mean(data[[name]])
    }
    if (mode == "sparse") {
      data[["y"]] = data[["wastage"]]
      explicit = y ~ log_outside_salary:factor(experience_group) +
        log_teacher_salary:factor(experience_group) + factor(spine_band_year) +
        factor(district) + factor(Experience) + factor(AGE) + factor(SEX) + factor(FTPTWK)
      full_df = 10L + sum(vapply(data[fes], function(x) length(unique(x)), integer(1))) - 5L
      covariance_df = full_df - length(unique(data[["district"]])) + 1L
      ref = teacher_sparse_reference(data, explicit, "district", covariance_df,
        initial_probability = data[["true_probability"]])
      stopifnot(ref[["model_df"]] == full_df, length(ref[["coef"]]) == 10L)
      saveRDS(ref, file.path(out, "sparse_reference.rds"))
      print(ref[c("nobs", "elapsed", "score", "slope_step", "model_df", "covariance_df")])
      return(invisible(NULL))
    }
    if (mode == "reference") {
      message("Starting independent strict fixest reference")
      elapsed = system.time({fit = dbglm_reference_feglm(formula, data,
        family = binomial("cloglog"), vcov = ~district, fixef.rm = "perfect_fit",
        etastart = binomial("cloglog")[["linkfun"]](data[["true_probability"]]), verbose = 1)})[["elapsed"]]
      ref = summarize_fixest(fit, data, elapsed)
      stopifnot(length(ref[["coef"]]) == 10L, ref[["converged"]])
      saveRDS(ref, file.path(out, "reference.rds"))
      print(ref[c("nobs", "elapsed", "refinements", "model_df", "covariance_df")])
      return(invisible(NULL))
    }
  }
  for (rep in seq_len(repetitions)) {
    gc()
    message("START ", mode, " repetition ", rep, " at ", Sys.time())
    if (mode == "fixest") {
      elapsed = system.time({fit = fixest::feglm(formula, data = data, family = binomial("cloglog"),
        vcov = ~district, fixef.rm = "perfect_fit", verbose = 1)})[["elapsed"]]
      result = summarize_fixest(fit, data, elapsed, rep)
    } else {
      elapsed = system.time({fit = dbglm(formula, conn = conn, table = table,
        family = binomial("cloglog"), vcov = ~district, row_id = "id",
        sample_table = retained, verbose = TRUE)})[["elapsed"]]
      keys = as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained, "ORDER BY id"))[["id"]])
      DBI::dbExecute(conn, paste(drop, retained))
      result = list(engine = mode, repetition = rep, elapsed = elapsed, coef = coef(fit),
        vcov = vcov(fit), logLik = as.numeric(logLik(fit)), nobs = nobs(fit),
        n_input = n, exits = as.numeric(summary[["exits"]]), keys = keys,
        model_df = fit[["n_params_model"]], covariance_df = fit[["n_params_vcov"]],
        iterations = fit[["iterations"]], converged = fit[["converged"]],
        score = fit[["score_max"]], projection = fit[["projection_residual"]],
        score_tol = fit[["control"]][["score_tol"]], fe_tol = fit[["control"]][["fe_tol"]])
      stopifnot(result[["score"]] <= result[["score_tol"]], result[["projection"]] <= result[["fe_tol"]])
    }
    stopifnot(length(result[["coef"]]) == 10L, result[["converged"]], result[["nobs"]] == n,
      identical(result[["keys"]], as.numeric(seq_len(n) - 1L)))
    saveRDS(result, file.path(out, paste0(mode, "_", rep, ".rds")))
    print(result[c("engine", "repetition", "elapsed", "iterations", "score", "projection")])
    fit = NULL
    result = NULL
    gc()
  }
}
run()
