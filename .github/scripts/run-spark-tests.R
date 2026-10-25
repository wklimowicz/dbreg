# Small real-Spark syntax/algorithm tier, deliberately not a speed benchmark.
args = commandArgs(TRUE)
out = if (length(args)) args[[1L]] else file.path(tempdir(), "dbreg-spark-validation")
home = Sys.getenv("SPARK_HOME", "/opt/spark")
release = file.path(home, "RELEASE")
if (!file.exists(release) || !grepl("^Spark 4\\.2\\.0 ", readLines(release, n = 1L))) {
  stop("This tier requires the final Spark 4.2.0 release at SPARK_HOME.", call. = FALSE)
}
if (!requireNamespace("sparklyr", quietly = TRUE)) stop("Install sparklyr for the local Spark tier.", call. = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out = normalizePath(out)
Sys.setenv(DBREG_TEST_SPARK = "true", DBREG_SPARK_TEST_OUT = out, SPARK_LOCAL_IP = "127.0.0.1")
pkgload::load_all()
fixest::setFixest_nthreads(2)
writeLines(c(readLines(release), capture.output(sessionInfo())), file.path(out, "environment.txt"))
sources = c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
  "inst/tinytest/helper_spark.R", "inst/tinytest/test_spark.R", ".github/scripts/run-spark-tests.R")
hashes = tools::md5sum(sources)
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)), file.path(out, "source-md5.csv"), row.names = FALSE)
started = Sys.time()
result = tryCatch(tinytest::run_test_file("inst/tinytest/test_spark.R"), error = identity)
if (inherits(result, "error")) {
  summary = data.frame(suite = "test_spark.R", assertions = 0L, failures = 1L, error = conditionMessage(result))
} else {
  saveRDS(result, file.path(out, "results.rds"))
  values = as.data.frame(result)
  summary = data.frame(suite = "test_spark.R", assertions = nrow(values),
    failures = if (nrow(values)) sum(is.na(values[["result"]]) | !values[["result"]]) else 1L,
    error = if (nrow(values)) "" else "No Spark assertions executed")
  print(result)
}
summary[["elapsed_seconds"]] = as.numeric(difftime(Sys.time(), started, units = "secs"))
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary)
if (summary[["failures"]] > 0L) stop("Local Spark validation failed; see summary.csv.", call. = FALSE)
