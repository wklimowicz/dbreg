# Production-code parity for the 88-slope model on the original synthetic rows.
# Run modes sequentially; preparation/reference fitting are not fit benchmarks.
# Rscript _benchmarks/teacher_cloglog_wide_validation.R INPUT OUT prepare|reference|duckdb|databricks [N]
# `prepare` regenerates the synthetic input on Databricks and saves it locally.
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("inst/tinytest/helper_glm_reference.R")
source("_benchmarks/teacher_cloglog_wide_compare.R")
args = commandArgs(TRUE)
input = args[[1]]
out = args[[2]]
mode = match.arg(args[[3]], c("prepare", "reference", "duckdb", "databricks"))
n = if (length(args) > 3L) as.integer(args[[4]]) else 5000000L
stopifnot(n >= 1000L, n <= 5000000L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
formula = wastage ~ log_outside_salary:factor(Experience) + log_teacher_salary:factor(Experience) |
  spine_band_year + district + Experience + AGE + SEX + FTPTWK
environment(formula) = baseenv()
# Shared DOUBLE arithmetic gives identical teacher wages on both engines.
teacher = teacher_cloglog_salary_sql()
reference_file = file.path(out, "reference.rds")
run = function() {
  fixest::setFixest_nthreads(2)
  remote = mode %in% c("prepare", "databricks")
  spill = file.path(out, "spill", paste0(mode, "-", Sys.getpid()))
  if (!remote) dir.create(spill, recursive = TRUE, showWarnings = FALSE)
  conn = if (remote) DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536) else
    DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB", temp_directory = spill)))
  on.exit(if (remote) DBI::dbDisconnect(conn) else DBI::dbDisconnect(conn, shutdown = TRUE))
  if (!remote) on.exit(unlink(spill, recursive = TRUE), add = TRUE)
  table = paste0("teacher_wide_", Sys.getpid())
  retained = paste0(table, "_retained")
  if (remote) {
    DBI::dbExecute(conn, "SET use_cached_result=false")
    query = paste0("SELECT *,", teacher, " FROM (", teacher_cloglog_sql(n), ") original")
  } else query = paste0("SELECT *,", teacher, " FROM read_parquet(",
    DBI::dbQuoteString(conn, file.path(input, "teachers.parquet")), ") WHERE id < ", n)
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", query))
  drop = if (remote) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
  on.exit(for (t in c(retained, table)) DBI::dbExecute(conn, paste(drop, t)), add = TRUE, after = FALSE)
  summary = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,SUM(wastage) AS exits,",
    "COUNT(DISTINCT Experience) AS experience_levels,COUNT(DISTINCT id) AS ids,MIN(id) AS lo,MAX(id) AS hi FROM", table))
  stopifnot(as.numeric(summary[["n"]]) == n, as.numeric(summary[["experience_levels"]]) == 44,
    as.numeric(summary[["ids"]]) == n, as.numeric(summary[["lo"]]) == 0, as.numeric(summary[["hi"]]) == n - 1L)
  metadata = c(capture.output(sessionInfo()), paste("mode:", mode, "rows:", n),
    "Synthetic only; DuckDB threads=4, fixest threads=2; default dbglm tolerances; result reuse disabled.",
    "DuckDB memory_limit=4GB; use an OUT directory on a disk with space for spilling.",
    "Elapsed times are validation durations, not isolated performance comparisons.",
    if (remote) capture.output(DBI::dbGetQuery(conn, "SELECT current_version() AS version")))
  writeLines(metadata, file.path(out, paste0(mode, "_environment.txt")))
  if (mode == "prepare") {
    dir.create(input, recursive = TRUE, showWarnings = FALSE)
    local = DBI::dbConnect(duckdb::duckdb(config = list(threads = "4", memory_limit = "4GB")))
    on.exit(DBI::dbDisconnect(local, shutdown = TRUE), add = TRUE)
    message("Downloading ", n, " synthetic rows for identical local/reference inputs")
    rows = DBI::dbGetQuery(conn, paste("SELECT * EXCEPT(log_teacher_salary) FROM", table))
    # Normalize ODBC integer64 columns before DuckDB's data-frame importer can
    # treat their stored bits as doubles. These synthetic IDs are below 2^53.
    for (name in names(rows)) if (inherits(rows[[name]], "integer64")) {
      rows[[name]] = as.numeric(rows[[name]])
    }
    stopifnot(identical(sort(rows[["id"]]), as.numeric(seq_len(n) - 1L)))
    DBI::dbWriteTable(local, "teachers", rows, temporary = TRUE)
    imported = DBI::dbGetQuery(local, "SELECT MIN(id) AS lo,MAX(id) AS hi,COUNT(DISTINCT id) AS n FROM teachers")
    stopifnot(imported[["lo"]] == 0, imported[["hi"]] == n - 1L, imported[["n"]] == n)
    DBI::dbExecute(local, paste0("COPY teachers TO ",
      DBI::dbQuoteString(local, file.path(input, "teachers.parquet")), " (FORMAT PARQUET)"))
    writeLines(teacher_cloglog_sql(n), file.path(input, "data_generation.sql"))
    write.csv(summary, file.path(out, "input_summary.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  if (mode == "reference") {
    d = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "ORDER BY id"))
    # Each wage is interacted with Experience, which is also an absorbed FE.
    # These independent constant shifts preserve every target slope and fitted
    # probability while conditioning fixest's own demeaning problem better.
    centers = vapply(d[c("log_outside_salary", "log_teacher_salary")], mean, numeric(1))
    for (name in names(centers)) d[[name]] = d[[name]] - centers[[name]]
    write.csv(data.frame(variable = names(centers), center = centers), file.path(out, "reference_centers.csv"), row.names = FALSE)
    # Known generating probabilities provide an independent starting point.
    # No dbglm result enters this reference; refinement/stability gates remain
    # identical to the other fixest comparisons. This is validation, not a
    # timed comparison against fixest's default initialization.
    initial_eta = binomial("cloglog")[["linkfun"]](d[["true_probability"]])
    message("Starting independent fixest reference: ", n, " rows, 88 slopes")
    elapsed = system.time({ref = dbglm_reference_feglm(formula, d, family = binomial("cloglog"),
      vcov = ~district, fixef.rm = "perfect_fit", fixef.tol = 1e-10,
      etastart = initial_eta, verbose = 1)})[["elapsed"]]
    stopifnot(length(coef(ref)) == 88L, isTRUE(ref[["convStatus"]]))
    keys = as.numeric(d[["id"]][fixest::obs(ref)])
    reference = list(coef = coef(ref), vcov = vcov(ref), logLik = as.numeric(logLik(ref)), nobs = nobs(ref),
      n_input = n, exits = as.numeric(summary[["exits"]]), keys = keys, model_df = ref[["nparams"]],
      covariance_df = attr(vcov(ref, attr = TRUE), "df.K"), elapsed = elapsed,
      refinements = attr(ref, "reference_refinements"), reference_engine = "fixest::feglm", fixef_tol = 1e-10)
    saveRDS(reference, reference_file)
    write.csv(data.frame(term = names(coef(ref)), coefficient = coef(ref), se = sqrt(diag(vcov(ref)))),
      file.path(out, "reference_coefficients.csv"), row.names = FALSE)
    message("Reference complete: ", elapsed, " seconds; nobs=", nobs(ref))
    return(invisible(NULL))
  }
  message("Starting ", mode, " wide fit")
  elapsed = system.time({fit = dbglm(formula, conn = conn, table = table, family = binomial("cloglog"),
    vcov = ~district, row_id = "id", sample_table = retained, verbose = TRUE)})[["elapsed"]]
  keys = as.numeric(DBI::dbGetQuery(conn, paste("SELECT id FROM", retained, "ORDER BY id"))[["id"]])
  bundle = list(fit = fit, keys = keys, n_input = n,
    exits = as.numeric(summary[["exits"]]), elapsed = elapsed, engine = mode)
  saveRDS(bundle, file.path(out, paste0(mode, "_pending_reference.rds")))
  if (file.exists(reference_file)) {
    teacher_wide_compare(bundle, readRDS(reference_file), out)
  } else {
    message("Fit saved. No fixest reference is available; fixest parity has NOT been checked.")
  }

}
run()
