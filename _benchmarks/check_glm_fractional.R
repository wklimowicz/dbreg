# Reproduce fractional-weight validation with durable results and source identity.
# Rscript _benchmarks/check_glm_fractional.R local|full|native OUTPUT_DIRECTORY [fractional|prediction|linear_prediction]
args = commandArgs(TRUE)
if (!length(args) %in% 2:3 || !args[[1L]] %in% c("local", "full", "native")) {
  stop("Usage: Rscript _benchmarks/check_glm_fractional.R local|full|native OUTPUT_DIRECTORY [fractional|prediction|linear_prediction]")
}
suite = if (length(args) == 3L) match.arg(args[[3L]], c("fractional", "prediction", "linear_prediction")) else "fractional"
suite_file = switch(suite, fractional = "test_glm_fractional", prediction = "test_prediction_state", linear_prediction = "test_linear_prediction")
mode = args[[1L]]
out = args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all()
fixest::setFixest_nthreads(2)
Sys.setenv(DBREG_TEST_DATABRICKS = if (mode == "native") "true" else "false")
writeLines(c(paste("UTC:", format(Sys.time(), tz = "UTC")),
  paste("Revision:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  paste("Fixest:", as.character(utils::packageVersion("fixest"))),
  paste("Native statement timeout:", getOption("dbreg.test_statement_timeout", 300L)),
  capture.output(sessionInfo())), file.path(out, paste0(mode, "-environment.txt")))
sources = c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
  list.files("inst/tinytest", full.names = TRUE),
  ".github/scripts/run-databricks-tests.R", "_benchmarks/check_glm_fractional.R")
hashes = tools::md5sum(sources)
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)),
  file.path(out, paste0(mode, "-source-md5.csv")), row.names = FALSE)
started = Sys.time()
result = tryCatch({
  if (mode == "full") tinytest::test_all(".") else tinytest::run_test_file(
    file.path("inst/tinytest", paste0(suite_file,
      if (mode == "native") "_live", ".R")))
}, error = identity)
elapsed = as.numeric(difftime(Sys.time(), started, units = "secs"))
if (inherits(result, "error")) {
  record = data.frame(mode = mode, assertions = 0L, failures = 1L,
    elapsed_seconds = elapsed, error = conditionMessage(result))
} else {
  saveRDS(result, file.path(out, paste0(mode, ".rds")))
  assertions = as.data.frame(result)
  record = data.frame(mode = mode, assertions = nrow(assertions),
    failures = if (nrow(assertions)) sum(is.na(assertions[["result"]]) | !assertions[["result"]]) else 1L,
    elapsed_seconds = elapsed, error = if (nrow(assertions)) "" else "Suite executed zero assertions")
  print(summary(result))
}
write.csv(record, file.path(out, paste0(mode, "-summary.csv")), row.names = FALSE)
print(record)
if (record[["failures"]] > 0L) stop("Validation failed; see saved artifacts.", call. = FALSE)
