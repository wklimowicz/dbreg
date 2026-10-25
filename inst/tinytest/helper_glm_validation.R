# Shared native/local validation plumbing. Only synthetic rows enter a backend.
dbglm_validation_connection = function(live = FALSE) {
  if (!live) return(DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "", threads = "2"))))
  timeout = getOption("dbreg.test_statement_timeout", 300L)
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) ||
    timeout < 1 || timeout > 172800 || timeout != floor(timeout)) {
    stop("Native test statement timeout must be an integer from 1 to 172800 seconds.", call. = FALSE)
  }
  conn = DBI::dbConnect(odbc::databricks(),
    driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
    httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE, DefaultStringColumnLength = 65536)
  ready = FALSE
  on.exit(if (!ready) DBI::dbDisconnect(conn))
  DBI::dbExecute(conn, "SET use_cached_result=false")
  # Bound native execution on the server, including calls that block in ODBC.
  # https://docs.databricks.com/aws/en/sql/language-manual/parameters/statement_timeout
  DBI::dbExecute(conn, paste("SET STATEMENT_TIMEOUT =", sprintf("%.0f", timeout)))
  ready = TRUE
  conn
}

dbglm_validation_tables = function(conn) {
  if (isTRUE(dbreg:::detect_backend(conn)[["is_databricks"]])) {
    tables = DBI::dbGetQuery(conn, "SHOW TABLES")
    names = tables[["tableName"]]
    benchmark = grepl("^dbreg_bench_[0-9a-f]{24}$", names)
    # Managed inputs from other benchmark connections are workspace-global.
    # Compare this session's temporary objects, and still catch accidental
    # persistent fit objects under every production-owned naming prefix.
    temporary = !is.na(tables[["isTemporary"]]) & as.logical(tables[["isTemporary"]])
    owned = grepl("^(dbglm_|dbreg_)", names) & !benchmark
    sort(names[temporary | owned])
  } else sort(DBI::dbListTables(conn))
}

dbglm_validation_drop = function(conn, table) {
  sql = if (isTRUE(dbreg:::detect_backend(conn)[["is_databricks"]])) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
  DBI::dbExecute(conn, paste(sql, DBI::dbQuoteIdentifier(conn, table)))
  invisible(NULL)
}

dbglm_validation_upload = function(conn, data) {
  table = basename(tempfile("glm_validation_"))
  if (isTRUE(dbreg:::detect_backend(conn)[["is_databricks"]])) {
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", dbreg:::sql_values_frame(conn, data)))
  } else DBI::dbWriteTable(conn, table, data, temporary = TRUE)
  table
}

dbglm_validation_compare = function(fit, reference, info = "", terms = names(coef(reference))) {
  reference_coef = coef(reference)[terms]
  reference_vcov = vcov(reference)[terms, terms, drop = FALSE]
  reference_se = sqrt(diag(reference_vcov))
  expect_identical(names(coef(fit)), names(reference_coef), info = info)
  expect_true(all(is.finite(coef(fit))) && all(is.finite(vcov(fit))), info = info)
  expect_true(all(abs(coef(fit) - reference_coef) <= 1e-7 + 1e-6 * abs(reference_coef)), info = info)
  expect_true(all(abs(vcov(fit) - reference_vcov) <= 1e-8 + 1e-5 * abs(reference_vcov)), info = info)
  expect_true(all(abs(sqrt(diag(vcov(fit))) - reference_se) <=
    1e-8 + 1e-6 * reference_se), info = info)
  expect_equal(nobs(fit), nobs(reference), info = info)
  expect_true(abs(as.numeric(logLik(fit) - logLik(reference))) / nobs(reference) < 1e-9, info = info)
  expect_equal(fit[["n_params_model"]], reference[["nparams"]], info = info)
  expect_equal(fit[["n_params_vcov"]], attr(vcov(reference, attr = TRUE), "df.K"), info = info)
  expect_true(fit[["converged"]] && fit[["score_max"]] <= fit[["control"]][["score_tol"]], info = info)
  expect_true(fit[["projection_residual"]] <= fit[["control"]][["fe_tol"]], info = info)
}
