# End-to-end comparison on identical synthetic data. Run from the package root:
# Rscript _benchmarks/backend_review.R [output-directory]
# Requires configured DATABRICKS_* variables and the Databricks ODBC driver.
# Only session-scoped remote objects are created; source tables are untouched.
pkgload::load_all()
args = commandArgs(trailingOnly = TRUE)
out = if (length(args)) args[[1]] else "."
dir.create(out, recursive = TRUE, showWarnings = FALSE)
local_conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "4"))
remote_conn = DBI::dbConnect(odbc::databricks(),
  driver = Sys.getenv("DATABRICKS_ODBC_DRIVER", "/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
  httpPath = Sys.getenv("DATABRICKS_SQL_PATH"), useNativeQuery = TRUE,
  DefaultStringColumnLength = 65536)
run = function() {
  on.exit(DBI::dbDisconnect(local_conn, shutdown = TRUE))
  on.exit(DBI::dbDisconnect(remote_conn), add = TRUE)
  options(dbreg.ap_tol = 1e-7)
  # Disable Databricks result reuse; warm warehouse startup is outside timers.
  DBI::dbExecute(remote_conn, "SET use_cached_result = false")
  metadata = list(session = sessionInfo(), duckdb_threads = 4,
    fixest_threads = fixest::getFixest_nthreads(),
    remote_version = DBI::dbGetQuery(remote_conn, "SELECT current_version() AS version"),
    linear_repetitions = 3L, glm_repetitions = 1L)
  saveRDS(metadata, file.path(out, "backend_review_metadata.rds"))
  state = new.env(parent = emptyenv())
  state[["rows"]] = list()
  record = function(n, workload, engine, repetition, thunk, reference) {
    gc()
    fit = NULL
    state[["error"]] = ""
    elapsed = system.time(tryCatch({fit = thunk()}, error = function(e) state[["error"]] = conditionMessage(e)))[["elapsed"]]
    state[["rows"]][[length(state[["rows"]]) + 1L]] = data.frame(n = n, workload = workload,
      engine = engine, repetition = repetition, elapsed = elapsed,
      coef_error = if (is.null(fit)) NA_real_ else max(abs(coef(fit) - coef(reference))),
      se_error = if (is.null(fit)) NA_real_ else max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(reference))))),
      iterations = if (length(fit[["iterations"]]) == 1L) fit[["iterations"]] else NA_integer_,
      method = if (!is.null(fit[["projection_method"]])) fit[["projection_method"]] else "",
      error = state[["error"]])
    results = do.call(rbind, state[["rows"]])
    write.csv(results, file.path(out, "backend_review.csv"), row.names = FALSE)
    print(utils::tail(results, 1), row.names = FALSE)
    flush.console()
  }
  for (n in c(100000L, 1000000L)) {
    table = paste0("dbreg_review_", Sys.getpid(), "_", n)
    # Generate remotely, then copy the resulting doubles to DuckDB and fixest.
    # This avoids unsupported ODBC parameter-array uploads in native mode.
    sql = paste0("CREATE TEMP TABLE ", table, " AS WITH raw AS (SELECT id, ",
      "CAST(pmod(id*17+floor(id/101),100)+1 AS INT) AS g, ",
      "CAST(pmod(id*13+floor(id/37),20)+1 AS INT) AS h, ",
      "CAST(pmod(id*7+floor(id/23),5)+1 AS INT) AS k, ",
      "CAST(pmod(id*11+floor(id/17),3)+1 AS INT) AS j, ",
      "(sin(id*12.9898)*43758.5453-floor(sin(id*12.9898)*43758.5453)-0.5)*4 AS x, ",
      "(sin(id*78.233)*19341.167-floor(sin(id*78.233)*19341.167)-0.5)*4 AS z, ",
      "sin(id*37.719)*9371.17-floor(sin(id*37.719)*9371.17) AS u, ",
      "sin(id*19.913)*7193.31-floor(sin(id*19.913)*7193.31) AS v ",
      "FROM range(", n, ")), pred AS (SELECT *, -.5+.4*x-.2*z+.3*sin(g)+.2*cos(h) AS eta FROM raw) ",
      "SELECT id,g,h,k,j,x,z,eta+(u-.5)*3.464101615 AS y, ",
      "CAST(v < 1-exp(-exp(eta)) AS INT) AS binary FROM pred")
    preparation = system.time({
      DBI::dbExecute(remote_conn, sql)
      d = DBI::dbGetQuery(remote_conn, paste("SELECT * FROM", table, "ORDER BY id"))
      DBI::dbWriteTable(local_conn, table, d, temporary = TRUE)
    })[["elapsed"]]
    cat("Prepared identical", n, "row inputs in", preparation, "seconds (excluded)\n")
    specs = list(
      moments_hc1 = list(fml = y ~ x + z, strategy = "moments", vcov = "hc1"),
      demean_cluster = list(fml = y ~ x + z | g + h, strategy = "demean", vcov = ~h),
      cloglog_4fe = list(fml = binary ~ x + z | g + h + k + j, family = binomial("cloglog"), vcov = ~h))
    for (name in names(specs)) {
      if (name == "cloglog_4fe" && n > 100000L) next
      spec = specs[[name]]
      is_glm = !is.null(spec[["family"]])
      reference_args = spec
      reference_args[["strategy"]] = NULL
      reference_args[["data"]] = d
      reference_args[["fixef.rm"]] = "none"
      if (is_glm) {
        reference_args[["glm.tol"]] = 1e-12
        reference_args[["glm.iter"]] = 100L
      }
      ref_fun = if (is_glm) fixest::feglm else fixest::feols
      reference = do.call(ref_fun, reference_args)
      for (rep in seq_len(if (is_glm) 1L else 3L)) {
        record(n, name, "fixest", rep, function() do.call(ref_fun, reference_args), reference)
        # Alternate backend order to reduce systematic run-order effects.
        engines = if (rep %% 2) c("DuckDB", "Databricks") else c("Databricks", "DuckDB")
        for (engine in engines) {
          conn = if (engine == "DuckDB") local_conn else remote_conn
          record(n, name, engine, rep,
            function() do.call(if (is_glm) dbglm else dbreg, c(spec, list(conn = conn, table = table))), reference)
        }
      }
    }
    DBI::dbRemoveTable(local_conn, table)
    DBI::dbExecute(remote_conn, paste("DROP TEMP TABLE", DBI::dbQuoteIdentifier(remote_conn, table)))
  }
}
run()
