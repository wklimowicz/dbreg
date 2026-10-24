# Rscript THIS DATABASE_FILE OUTPUT_PREFIX
# Reuses the SQL-only paired-noise fixture prepared by
# ../two_fe_prediction_20261024/large-state.R; does not collect its panel.
args = commandArgs(TRUE)
stopifnot(length(args) == 2L, file.exists(args[[1L]]))
pkgload::load_all()
options(dbreg.prediction_max_entries = 2e6)
conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[1L]])
initial = sort(DBI::dbListTables(conn))
sources = c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE))
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))), paste0(args[[2L]], "-source-md5.csv"), row.names = FALSE)
writeLines(c(paste("UTC:", format(Sys.time(), tz = "UTC")), capture.output(sessionInfo())), paste0(args[[2L]], "-environment.txt"))
fit = dbreg(y ~ x + z | g + h, conn = conn, table = "panel", strategy = "demean", retain = "prediction")
stopifnot(identical(fit[["fitted_state"]][["prediction"]][["storage"]], "database"))
fitted_tables = sort(DBI::dbListTables(conn))
sql_prediction_metrics = new.env(parent = emptyenv())
sql_prediction_metrics[["active"]] = FALSE
trace("dbGetQuery", where = asNamespace("DBI"), print = FALSE, exit = quote({
  if (isTRUE(sql_prediction_metrics[["active"]])) {
    value = returnValue()
    sql_prediction_metrics[["rows"]] = max(sql_prediction_metrics[["rows"]], nrow(value))
    sql_prediction_metrics[["bytes"]] = sql_prediction_metrics[["bytes"]] + as.numeric(object.size(value))
  }
}))
records = list()
for (i in 0:3) {
  gc()
  sql_prediction_metrics[["rows"]] = 0L
  sql_prediction_metrics[["bytes"]] = 0
  sql_prediction_metrics[["active"]] = TRUE
  elapsed = system.time({output = predict(fit, conn = conn, table = "panel", row_id = "id")})[["elapsed"]]
  sql_prediction_metrics[["active"]] = FALSE
  maximum = sql_prediction_metrics[["rows"]]
  bytes = sql_prediction_metrics[["bytes"]]
  # Correctness stays in SQL; keyed join independently checks every prediction.
  quoted = as.character(DBI::dbQuoteIdentifier(conn, output[["table"]]))
  check = DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n,MAX(ABS(o.prediction-(2+.3*p.x-.2*p.z+SIN(p.g*.05)+COS(p.h*.03)))) AS error FROM ", quoted, " o JOIN panel p ON o.id=p.id"))
  stopifnot(as.numeric(check[["n"]]) == 1e6, check[["error"]] < 1e-7,
    maximum <= 1L, length(setdiff(DBI::dbListTables(conn), fitted_tables)) == 1L)
  stopifnot(nrow(head(output, 6L)) == 6L)
  if (i > 0L) records[[i]] = data.frame(repetition = i, rows = as.numeric(check[["n"]]),
    elapsed_seconds = elapsed, maximum_collected_rows = maximum, collected_object_bytes = bytes,
    max_prediction_error = check[["error"]], result_object_bytes = as.numeric(object.size(output)))
  close(output)
  stopifnot(identical(sort(DBI::dbListTables(conn)), fitted_tables))
  message("Million-row keyed output ", i, ": ", elapsed, "s")
}
untrace("dbGetQuery", where = asNamespace("DBI"))
close(fit)
stopifnot(identical(sort(DBI::dbListTables(conn)), initial), DBI::dbIsValid(conn))
DBI::dbDisconnect(conn, shutdown = TRUE)
write.csv(do.call(rbind, records), paste0(args[[2L]], ".csv"), row.names = FALSE)
