# Rscript THIS prepare DBFILE
# Rscript THIS run DBFILE none|prediction OUTPUT_PREFIX
args = commandArgs(TRUE)
stopifnot(length(args) %in% c(2L, 4L))
if (args[[1L]] == "prepare") {
  stopifnot(!file.exists(args[[2L]]))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[2L]])
  DBI::dbExecute(conn, paste0("CREATE TABLE panel AS WITH indices AS (SELECT i AS id,CAST(FLOOR((i-1)/2.0) AS BIGINT) AS j FROM range(1,1000001) t(i)),",
    "base AS (SELECT id,MOD(j,50000) AS g,CAST(MOD(hash(j*37),50000) AS BIGINT) AS h,",
    "SIN(j*.71) AS x,COS(j*.31) AS z FROM indices) ",
    "SELECT *,2+.3*x-.2*z+SIN(g*.05)+COS(h*.03)+CASE WHEN MOD(id,2)=0 THEN .1 ELSE -.1 END AS y FROM base"))
  print(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n,COUNT(DISTINCT g) AS groups_g,COUNT(DISTINCT h) AS groups_h FROM panel"))
  DBI::dbDisconnect(conn, shutdown = TRUE)
} else {
  stopifnot(args[[1L]] == "run", args[[3L]] %in% c("none", "prediction"))
  pkgload::load_all()
  mode = args[[3L]]
  out = args[[4L]]
  options(dbreg.prediction_max_entries = 2e6)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[2L]])
  metrics = new.env(parent = emptyenv())
  metrics[["maximum_rows"]] = 0L
  metrics[["total_bytes"]] = 0
  original_query = get("dbglm_get_query", asNamespace("dbreg"))
  measured_query = function(...) {
    value = original_query(...)
    if (identical(metrics[["phase"]], "recovery")) {
      metrics[["maximum_rows"]] = max(metrics[["maximum_rows"]], nrow(value))
      metrics[["total_bytes"]] = metrics[["total_bytes"]] + as.numeric(object.size(value))
    }
    value
  }
  assignInNamespace("dbglm_get_query", measured_query, "dbreg")
  original_recovery = get("build_database_two_fe_state", asNamespace("dbreg"))
  measured_recovery = function(...) {
    metrics[["phase"]] = "recovery"
    on.exit({metrics[["phase"]] = "fit"})
    original_recovery(...)
  }
  assignInNamespace("build_database_two_fe_state", measured_recovery, "dbreg")
  sources = list.files("R", full.names = TRUE)
  write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))), paste0(out, "-source-md5.csv"), row.names = FALSE)
  writeLines(c(paste("Mode:", mode), "N=1,000,000; 50,000 levels per FE block (approximately); two slopes; default prediction budget=2e6.",
    "DuckDB threads=2, memory_limit=2GB; process RSS measured externally.", capture.output(sessionInfo())), paste0(out, "-environment.txt"))
  initial_tables = sort(DBI::dbListTables(conn))
  records = list()
  values = list()
  for (i in 0:3) {
    metrics[["maximum_rows"]] = 0L
    metrics[["total_bytes"]] = 0
    gc()
    elapsed = system.time({fit = dbreg(y ~ x + z | g + h, conn = conn,
      table = "panel", strategy = "demean", retain = mode)})[["elapsed"]]
    fit_maximum_rows = metrics[["maximum_rows"]]
    fit_total_bytes = metrics[["total_bytes"]]
    if (mode == "prediction") {
      stopifnot(identical(fit[["fitted_state"]][["resources"]][["storage"]], "database"))
      data = DBI::dbGetQuery(conn, "SELECT x,z,g,h FROM panel ORDER BY id LIMIT 100")
      prediction = predict(fit, data)
      # Paired rows have identical regressors/FEs and opposite noise, so the
      # noise is exactly orthogonal to every model column. The true fitted
      # predictor is known without collecting training observations into R.
      reference = 2 + .3*data[["x"]] - .2*data[["z"]] + sin(data[["g"]]*.05) + cos(data[["h"]]*.03)
      stopifnot(all(is.finite(prediction)), max(abs(prediction - reference)) < 1e-7)
      dictionary = fit[["fitted_state"]][["specification"]][["dictionaries"]][[1L]][["table"]]
      stopifnot(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", dictionary))[["n"]]) == 50000)
      stopifnot(length(fit[["fitted_state"]][["specification"]][["dictionaries"]]) == 2L,
        as.numeric(object.size(fit[["fitted_state"]])) < 128*1024)
    }
    if (i > 0L) {
      records[[i]] = data.frame(repetition = i, mode = mode, elapsed_seconds = elapsed,
        nobs = nobs(fit), state_bytes = as.numeric(object.size(fit[["fitted_state"]])),
        maximum_prediction_fit_query_rows = fit_maximum_rows, prediction_fit_query_bytes = fit_total_bytes)
      values[[i]] = list(coef = coef(fit), vcov = vcov(fit), nobs = nobs(fit))
    }
    close(fit)
    stopifnot(identical(sort(DBI::dbListTables(conn)), initial_tables), DBI::dbIsValid(conn))
    message(mode, " repetition ", i, ": ", elapsed, "s")
  }
  write.csv(do.call(rbind, records), paste0(out, ".csv"), row.names = FALSE)
  saveRDS(values, paste0(out, "-values.rds"))
  DBI::dbDisconnect(conn, shutdown = TRUE)
}
