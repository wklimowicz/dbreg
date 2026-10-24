# Called from the repository root by CI or `make test-databricks`.
# Missing credentials/driver or an empty suite fails; it never counts as a pass.
# Optional arguments: output directory, then one or more known suite filenames.
args = commandArgs(TRUE)
out = if (length(args)) args[[1]] else file.path(tempdir(), "dbreg-databricks-validation")
required = c("DATABRICKS_HOST", "DATABRICKS_TOKEN", "DATABRICKS_SQL_PATH")
missing = required[!nzchar(Sys.getenv(required))]
if (length(missing)) stop("Missing Databricks configuration: ", paste(missing, collapse = ", "), call. = FALSE)
driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so")
if (!file.exists(driver)) stop("Databricks ODBC driver is unavailable: ", driver, call. = FALSE)
Sys.setenv(DBREG_TEST_DATABRICKS = "true")
pkgload::load_all()
fixest::setFixest_nthreads(2)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("UTC:", format(Sys.time(), tz = "UTC")), capture.output(sessionInfo())), file.path(out, "environment.txt"))
writeLines(c(paste("dbreg revision:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  paste("fixest version:", as.character(utils::packageVersion("fixest"))),
  paste("native statement timeout:", getOption("dbreg.test_statement_timeout", 300L)),
  "Working-tree contents are identified by source-md5.csv; revision alone may include uncommitted changes."),
  file.path(out, "versions.txt"))
sources = c("DESCRIPTION", "NAMESPACE", ".github/scripts/run-databricks-tests.R",
  list.files("R", full.names = TRUE), list.files("inst/tinytest", full.names = TRUE))
hashes = tools::md5sum(sources)
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)), file.path(out, "source-md5.csv"), row.names = FALSE)
suites = c("test_glm_databricks_live.R", "test_glm_databricks_sample.R", "test_fe_edges_live.R",
  "test_ap_moments_live.R", "test_glm_stress_live.R", "test_glm_separation_live.R", "test_glm_wide_live.R", "test_linear_design_live.R",
  "test_glm_fractional_live.R", "test_sample_identity_live.R", "test_prediction_state_live.R", "test_prediction_snapshot_live.R",
  "test_linear_prediction_live.R", "test_database_prediction_live.R")
if (length(args) > 1L) {
  requested = args[-1L]
  if (anyDuplicated(requested) || any(!requested %in% suites)) {
    stop("Specify unique known native suite filenames after the output directory.", call. = FALSE)
  }
  suites = requested
}
records = list()
for (file in suites) {
  message("Starting native suite: ", file)
  started = Sys.time()
  result = tryCatch(tinytest::run_test_file(file.path("inst/tinytest", file)), error = identity)
  if (inherits(result, "error")) {
    record = data.frame(suite = file, assertions = 0L, failures = 1L, error = conditionMessage(result))
  } else {
    saveRDS(result, file.path(out, paste0(file, ".rds")))
    assertions = as.data.frame(result)
    failures = sum(is.na(assertions[["result"]]) | !assertions[["result"]])
    record = data.frame(suite = file, assertions = nrow(assertions),
      failures = if (nrow(assertions)) failures else 1L,
      error = if (nrow(assertions)) "" else "Native suite executed zero assertions")
    print(result)
  }
  record[["elapsed_seconds"]] = as.numeric(difftime(Sys.time(), started, units = "secs"))
  records[[length(records) + 1L]] = record
  write.csv(do.call(rbind, records), file.path(out, "summary.csv"), row.names = FALSE)
}
summary = do.call(rbind, records)
print(summary)
if (any(summary[["failures"]] > 0L)) stop("Native Databricks validation failed; see summary.csv.", call. = FALSE)
