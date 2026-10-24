# Same-data, same-runtime comparison of the pre-compatibility and current code.
# prepare DBFILE N
# run REPOSITORY DBFILE default|integer|row_id|sample|fraction|prediction OUTPUT_PREFIX REPEATS LABEL
args = commandArgs(TRUE)
if (length(args) < 1L) stop("Expected prepare or run.")
if (args[[1L]] == "prepare") {
  if (length(args) != 3L) stop("prepare DBFILE N")
  n = as.numeric(args[[3L]])
  stopifnot(is.finite(n), n > 0, n == floor(n), !file.exists(args[[2L]]))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[2L]])
  DBI::dbExecute(conn, paste0("CREATE TABLE perf_input AS WITH base AS (SELECT i::BIGINT AS id,",
    "SIN(i*.001) AS x,COS(i*.00037) AS z,MOD(i-1,1000) AS g,",
    "CAST(MOD(FLOOR((i-1)/1000.0),50) AS BIGINT) AS h FROM range(1,", sprintf("%.0f", n + 1), ") t(i)) ",
    "SELECT *, .5+.3*x-.2*z+.01*g+.03*h+SIN(id*.17) AS y_lin,",
    "CASE WHEN CAST(MOD(hash(id),10000) AS DOUBLE)<10000*(1-EXP(-EXP(-.3+.2*x-.15*z+.1*SIN(g)+.1*COS(h)))) ",
    "THEN 1.0 ELSE 0.0 END AS y_bin,CAST(1+MOD(id,3) AS DOUBLE) AS w_int,",
    "CAST(.25*(1+MOD(id,3)) AS DOUBLE) AS w_frac FROM base"))
  print(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM perf_input"))
  DBI::dbDisconnect(conn, shutdown = TRUE)
} else if (args[[1L]] == "run") {
  if (length(args) != 7L) stop("run REPOSITORY DBFILE MODE OUTPUT_PREFIX REPEATS LABEL")
  repo = normalizePath(args[[2L]])
  mode = match.arg(args[[4L]], c("default", "integer", "row_id", "sample", "fraction", "prediction"))
  out = args[[5L]]
  repeats = as.integer(args[[6L]])
  # The final argument is a descriptive revision label, separate from HEAD.
  label = args[[7L]]
  stopifnot(repeats >= 1L)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  pkgload::load_all(repo)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")),
    dbdir = args[[3L]], bigint = "integer64")
  metrics = new.env(parent = emptyenv())
  metrics[["ll_seconds"]] = 0
  metrics[["ll_calls"]] = 0L
  ns = asNamespace("dbreg")
  if (exists("dbglm_binomial_loglik", ns, inherits = FALSE)) {
    original_ll = get("dbglm_binomial_loglik", ns)
    timed_ll = function(...) {
      started = proc.time()[["elapsed"]]
      value = original_ll(...)
      metrics[["ll_seconds"]] = metrics[["ll_seconds"]] + proc.time()[["elapsed"]] - started
      metrics[["ll_calls"]] = metrics[["ll_calls"]] + 1L
      value
    }
    assignInNamespace("dbglm_binomial_loglik", timed_ll, "dbreg")
  } else if (mode == "fraction") stop("Fractional weights are unsupported in this revision.")
  revision = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
  environment = c(paste("Label:", label), paste("Revision:", revision), paste("Mode:", mode),
    "DuckDB threads=2; memory_limit=2GB; persisted synthetic input; no panel collected into R.",
    capture.output(DBI::dbGetQuery(conn, "SELECT current_setting('threads') AS threads,current_setting('memory_limit') AS memory_limit")),
    "Only the fractional likelihood reporter is timed in-process; numerical code is unchanged.",
    "Two warm-up fits per model are excluded; measured repetitions follow.",
    capture.output(sessionInfo()))
  writeLines(trimws(environment, which = "right"), paste0(out, "-environment.txt"))
  files = list.files(file.path(repo, "R"), full.names = TRUE)
  hashes = tools::md5sum(files)
  write.csv(data.frame(file = substring(names(hashes), nchar(repo) + 2L), md5 = unname(hashes)),
    paste0(out, "-source-md5.csv"), row.names = FALSE)
  models = c("ols", "ols_fe", "glm", "glm_fe")
  if (mode %in% c("fraction", "prediction")) models = c("glm", "glm_fe")
  rows = list()
  run_one = function(model, repetition) {
    glm = startsWith(model, "glm")
    fe = endsWith(model, "_fe")
    fml = stats::as.formula(paste(if (glm) "y_bin" else "y_lin", "~x+z", if (fe) "|g+h" else ""))
    fit_args = list(fml = fml, conn = conn, table = "perf_input")
    if (glm) fit_args[["family"]] = stats::binomial("cloglog")
    if (mode %in% c("row_id", "sample")) fit_args[["row_id"]] = "id"
    if (mode == "sample") fit_args[["sample_table"]] = "bench_retained"
    if (mode == "integer") fit_args[["weights"]] = "w_int"
    if (mode == "fraction") fit_args[["weights"]] = "w_frac"
    if (mode == "prediction") fit_args[["retain"]] = "prediction"
    metrics[["ll_seconds"]] = 0
    metrics[["ll_calls"]] = 0L
    gc()
    elapsed = system.time({fit = do.call(if (glm) dbglm else dbreg, fit_args)})[["elapsed"]]
    if (glm && !isTRUE(fit[["converged"]])) stop("Only completed converged fits are benchmarked.")
    if (mode == "prediction" && fe) {
      newdata = DBI::dbGetQuery(conn, "SELECT x,z,g,h FROM perf_input ORDER BY id LIMIT 100")
      prediction = stats::predict(fit, newdata)
      stopifnot(all(is.finite(prediction)),
        identical(stats::predict(unserialize(serialize(fit, NULL)), newdata), prediction))
    }
    record = data.frame(label = label, revision = revision, mode = mode, model = model,
      repetition = repetition, elapsed_seconds = elapsed, nobs = nobs(fit),
      strategy = fit[["strategy"]], iterations = if (glm) fit[["iterations"]] else NA_integer_,
      likelihood_seconds = metrics[["ll_seconds"]], likelihood_calls = metrics[["ll_calls"]],
      prediction_state_bytes = as.numeric(object.size(fit[["fitted_state"]])),
      prediction_recovery_error = if (is.null(fit[["fitted_state"]])) NA_real_ else fit[["fitted_state"]][["prediction"]][["recovery_error"]],
      weighted_loglik = if (glm) as.numeric(if (mode == "fraction") logLik(fit, type = "weighted") else logLik(fit)) else NA_real_,
      coef_names = as.character(jsonlite::toJSON(names(coef(fit)), auto_unbox = FALSE)),
      coef_values = as.character(jsonlite::toJSON(unname(coef(fit)), digits = 17, auto_unbox = FALSE)),
      vcov_values = as.character(jsonlite::toJSON(as.numeric(vcov(fit)), digits = 17, auto_unbox = FALSE)))
    if (mode == "sample") {
      count = as.numeric(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM bench_retained")[["n"]])
      stopifnot(count == nobs(fit))
      DBI::dbRemoveTable(conn, "bench_retained")
    }
    record
  }
  for (model in models) {
    for (warmup in 1:2) invisible(run_one(model, 0L))
    for (repetition in seq_len(repeats)) {
      record = run_one(model, repetition)
      rows[[length(rows) + 1L]] = record
      write.csv(do.call(rbind, rows), paste0(out, ".csv"), row.names = FALSE)
      message(label, " ", mode, " ", model, " ", repetition, ": ", record[["elapsed_seconds"]], "s")
    }
  }
  DBI::dbDisconnect(conn, shutdown = TRUE)
} else stop("Expected prepare or run.")
