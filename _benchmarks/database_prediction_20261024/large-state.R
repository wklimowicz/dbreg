# Rscript THIS prepare DBFILE
# Rscript THIS run DBFILE none|prediction OUTPUT_PREFIX
args = commandArgs(TRUE)
stopifnot(length(args) %in% c(2L, 4L))
if (args[[1L]] == "prepare") {
  stopifnot(!file.exists(args[[2L]]))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[2L]])
  DBI::dbExecute(conn, paste0("CREATE TABLE panel AS WITH base AS (SELECT i AS id,MOD(i-1,400000) AS g,",
    "SIN(i*.71) AS x,COS(i*.31) AS z,SIN(i*.13) AS u,COS(i*.07) AS v FROM range(1,1000001) t(i)) ",
    "SELECT *,2+.3*x-.2*z+.1*u-.4*v+SIN(g*.05)+.1*SIN(id*.17) AS y FROM base"))
  print(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n,COUNT(DISTINCT g) AS groups FROM panel"))
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
    metrics[["maximum_rows"]] = max(metrics[["maximum_rows"]], nrow(value))
    metrics[["total_bytes"]] = metrics[["total_bytes"]] + as.numeric(object.size(value))
    value
  }
  assignInNamespace("dbglm_get_query", measured_query, "dbreg")
  sources = list.files("R", full.names = TRUE)
  write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))), paste0(out, "-source-md5.csv"), row.names = FALSE)
  writeLines(c(paste("Mode:", mode), "N=1,000,000; G=400,000; four slopes; default prediction budget=2e6.",
    "DuckDB threads=2, memory_limit=2GB; process RSS measured externally.", capture.output(sessionInfo())), paste0(out, "-environment.txt"))
  initial_tables = sort(DBI::dbListTables(conn))
  records = list()
  values = list()
  for (i in 0:3) {
    metrics[["maximum_rows"]] = 0L
    metrics[["total_bytes"]] = 0
    gc()
    elapsed = system.time({fit = dbreg(y ~ x + z + u + v | g, conn = conn,
      table = "panel", strategy = "demean", retain = mode)})[["elapsed"]]
    fit_maximum_rows = metrics[["maximum_rows"]]
    fit_total_bytes = metrics[["total_bytes"]]
    if (mode == "prediction") {
      stopifnot(identical(fit[["fitted_state"]][["resources"]][["storage"]], "database"))
      data = DBI::dbGetQuery(conn, "SELECT x,z,u,v,g FROM panel ORDER BY id LIMIT 100")
      prediction = predict(fit, data)
      beta = coef(fit)
      xb = paste(sprintf("(%.17g)*%s", beta, names(beta)), collapse = "+")
      reference = DBI::dbGetQuery(conn, paste0("SELECT ", xb, "+AVG(y-(", xb,
        ")) OVER (PARTITION BY g) AS prediction FROM panel ORDER BY id LIMIT 100"))[["prediction"]]
      stopifnot(all(is.finite(prediction)), max(abs(prediction - reference)) < 1e-7)
      dictionary = fit[["fitted_state"]][["specification"]][["dictionaries"]][[1L]][["table"]]
      stopifnot(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", dictionary))[["n"]]) == 400000)
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
