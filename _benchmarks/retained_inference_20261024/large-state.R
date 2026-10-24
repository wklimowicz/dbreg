# Rscript THIS DATABASE_FILE prediction|all OUTPUT_PREFIX
# Uses the preceding two_fe_prediction_20261024 paired-noise SQL fixture.
args = commandArgs(TRUE)
stopifnot(length(args) == 3L, file.exists(args[[1L]]), args[[2L]] %in% c("prediction", "all"))
pkgload::load_all()
mode = args[[2L]]
out = args[[3L]]
options(dbreg.prediction_max_entries = 2e6)
conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2", memory_limit = "2GB")), dbdir = args[[1L]])
original_tables = sort(DBI::dbListTables(conn))
sources = c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE))
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))), paste0(out, "-source-md5.csv"), row.names = FALSE)
writeLines(c(paste("UTC:", format(Sys.time(), tz = "UTC")), paste("Mode:", mode),
  capture.output(sessionInfo())), paste0(out, "-environment.txt"))
retained_metrics = new.env(parent = emptyenv())
retained_metrics[["active"]] = FALSE
trace("dbGetQuery", where = asNamespace("DBI"), print = FALSE, exit = quote({
  if (isTRUE(retained_metrics[["active"]])) {
    value = returnValue()
    retained_metrics[["rows"]] = max(retained_metrics[["rows"]], nrow(value))
    retained_metrics[["bytes"]] = retained_metrics[["bytes"]] + as.numeric(object.size(value))
  }
}))
original_capture = get("retain_score_sample", asNamespace("dbreg"))
measured_capture = function(...) {
  retained_metrics[["active"]] = TRUE
  on.exit({retained_metrics[["active"]] = FALSE})
  original_capture(...)
}
assignInNamespace("retain_score_sample", measured_capture, "dbreg")
records = list()
values = list()
for (i in 0:3) {
  gc()
  retained_metrics[["rows"]] = 0L
  retained_metrics[["bytes"]] = 0
  fit_time = system.time({fit = dbreg(y ~ x + z | g + h, conn = conn, table = "panel",
    strategy = "demean", row_id = "id", retain = mode)})[["elapsed"]]
  capture_rows = retained_metrics[["rows"]]
  capture_bytes = retained_metrics[["bytes"]]
  post_time = NA_real_
  post_rows = post_bytes = 0
  if (mode == "all") {
    state = fit[["fitted_state"]][["inference"]]
    stopifnot(state[["nobs"]] == 1e6, state[["nobs_input"]] == 1e6, capture_rows <= 1L,
      length(setdiff(DBI::dbListTables(conn), original_tables)) == 3L)
    retained_metrics[["rows"]] = 0L
    retained_metrics[["bytes"]] = 0
    retained_metrics[["active"]] = TRUE
    post_time = system.time({HC = vcov(fit, vcov = "hc1"); clustered = vcov(fit, vcov = ~g+h, ssc = list(K.fixef = "nonnested"))})[["elapsed"]]
    retained_metrics[["active"]] = FALSE
    post_rows = retained_metrics[["rows"]]
    post_bytes = retained_metrics[["bytes"]]
    # Every residual is +/-0.1 in this paired-noise fixture. HC1 and full-K
    # IID covariance therefore agree, without collecting X or residual rows.
    stopifnot(post_rows <= 1L, max(abs(HC-vcov(fit))) < 1e-12, all(is.finite(clustered)))
    output = predict(fit)
    table = as.character(DBI::dbQuoteIdentifier(conn, output[["table"]]))
    check = DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n,MAX(ABS(o.prediction-(2+.3*p.x-.2*p.z+SIN(p.g*.05)+COS(p.h*.03)))) AS error FROM ",
      table, " o JOIN panel p ON o.id=p.id"))
    stopifnot(as.numeric(check[["n"]]) == 1e6, check[["error"]] < 1e-7, nrow(head(output)) == 6L)
  }
  if (i > 0L) {
    records[[i]] = data.frame(repetition = i, mode = mode, elapsed_fit_seconds = fit_time,
      elapsed_postfit_seconds = post_time, capture_max_rows = capture_rows, capture_object_bytes = capture_bytes,
      postfit_max_rows = post_rows, postfit_object_bytes = post_bytes,
      model_object_bytes = as.numeric(object.size(fit)))
    values[[i]] = list(coef = coef(fit), vcov = vcov(fit), nobs = nobs(fit))
  }
  close(fit)
  if (mode == "all") {
    stopifnot(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]]) == 1e6)
    close(output)
  }
  stopifnot(identical(sort(DBI::dbListTables(conn)), original_tables), DBI::dbIsValid(conn))
  message(mode, " repetition ", i, ": ", fit_time, "s fit; ", post_time, "s post-fit")
}
untrace("dbGetQuery", where = asNamespace("DBI"))
DBI::dbDisconnect(conn, shutdown = TRUE)
write.csv(do.call(rbind, records), paste0(out, ".csv"), row.names = FALSE)
saveRDS(values, paste0(out, "-values.rds"))
