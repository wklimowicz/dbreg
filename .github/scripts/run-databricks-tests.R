# Called from the repository root by CI or `make test-databricks`.
# Missing credentials/driver or an empty suite fails; it never counts as a pass.
# Optional arguments: output directory, then known suite filenames or --group NAME.
# --list-suites validates and lists the inventory without credentials or a driver.
local({
  # A single top-level expression is parsed before any long-running test starts.
  # Source hashes below still reject edits, without incremental Rscript parsing.
  args = commandArgs(TRUE)
  out = if (length(args)) args[[1]] else file.path(tempdir(), "dbreg-databricks-validation")
  suites = readLines("inst/tinytest/native-suites.txt")
  if (!length(suites) || anyDuplicated(suites) || any(!grepl("^test_[a-z_]+[.]R$", suites)) ||
      any(!file.exists(file.path("inst/tinytest", suites)))) {
    stop("Invalid native suite inventory: every entry must name a unique existing test file.", call. = FALSE)
  }
  if (identical(args, "--list-suites")) {
    writeLines(suites)
    return(invisible(NULL))
  }

  groups = list(
    `glm-core` = c("test_glm_databricks_live.R", "test_glm_databricks_sample.R", "test_glm_fractional_live.R", "test_gaussian_live.R"),
    `glm-stress` = c("test_fe_edges_live.R", "test_glm_stress_live.R", "test_glm_separation_live.R", "test_glm_wide_live.R"),
    prediction = suites[grepl("prediction", suites)])
  groups[["linear-inference"]] = setdiff(suites, unlist(groups, use.names = FALSE))
  stopifnot(!anyDuplicated(unlist(groups)), setequal(unlist(groups), suites))
  groups[["smoke"]] = c("test_projection_transport_live.R", "test_glm_smoke_live.R", "test_ap_moments_live.R", "test_sql_weights_live.R")
  if (length(args) > 1L) {
    requested = args[-1L]
    if (length(requested) == 2L && requested[[1L]] == "--group") {
      if (!requested[[2L]] %in% names(groups)) stop("Unknown native suite group.", call. = FALSE)
      requested = groups[[requested[[2L]]]]
    }
    if (anyDuplicated(requested) || any(!requested %in% suites)) {
      stop("Specify unique known native suite filenames or --group NAME after the output directory.", call. = FALSE)
    }
    suites = requested
  }
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
  records = list()
  core_sources = c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE))
  loaded_core = tools::md5sum(core_sources)
  for (file in suites) {
    if (!identical(tools::md5sum(core_sources), loaded_core)) stop("Loaded package source changed during native validation; rerun from an immutable checkout.", call. = FALSE)
    suite_hashes = tools::md5sum(sources)
    write.csv(data.frame(file = names(suite_hashes), md5 = unname(suite_hashes)),
      file.path(out, paste0(file, ".source-md5.csv")), row.names = FALSE)
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
    if (!identical(tools::md5sum(sources), suite_hashes)) {
      record[["failures"]] = record[["failures"]]+1L
      record[["error"]] = paste(record[["error"]], "Validation source changed within this suite; provenance is unstable")
    }
    record[["elapsed_seconds"]] = as.numeric(difftime(Sys.time(), started, units = "secs"))
    records[[length(records) + 1L]] = record
    write.csv(do.call(rbind, records), file.path(out, "summary.csv"), row.names = FALSE)
  }
  summary = do.call(rbind, records)
  print(summary)
  if (any(summary[["failures"]] > 0L)) stop("Native Databricks validation failed; see summary.csv.", call. = FALSE)
})
