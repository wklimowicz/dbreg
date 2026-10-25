# Run each mode in a fresh /usr/bin/time process to include embedded DuckDB RSS.
args = commandArgs(TRUE)
stopifnot(length(args)==2L,args[[2L]] %in% c("baseline","lean","full"))
pkgload::load_all()
mode = args[[2L]]
root = paste0("_benchmarks/formula_i_20261024/large-",mode)
sources = c("DESCRIPTION","NAMESPACE",list.files("R",full.names=TRUE),"_benchmarks/formula_i_20261024/large-fit.R")
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),paste0(root,"-source-md5.csv"),row.names=FALSE)
writeLines(c(paste("UTC:",format(Sys.time(),tz="UTC")),capture.output(sessionInfo())),paste0(root,"-environment.txt"))
conn = DBI::dbConnect(duckdb::duckdb(config=list(threads="2",memory_limit="2GB")),dbdir=args[[1L]])
DBI::dbExecute(conn,"CREATE TEMP VIEW dbreg_i_perf AS SELECT *,CASE WHEN id%3=0 THEN 'a' WHEN id%3=1 THEN 'b' ELSE 'c' END AS f FROM panel")
original = sort(DBI::dbListTables(conn))
reference = dbreg(y~x+z+f,conn=conn,table="dbreg_i_perf",strategy="moments",row_id="id")
reference_ll = as.numeric(logLik(reference))
metrics = new.env(parent=emptyenv())
metrics[["active"]] = FALSE
trace("dbGetQuery",where=asNamespace("DBI"),print=FALSE,exit=quote({
 if(isTRUE(metrics[["active"]])) {
  metrics[["reads"]] = metrics[["reads"]]+1L
  value = returnValue()
  metrics[["rows"]] = max(metrics[["rows"]],nrow(value))
  metrics[["bytes"]] = metrics[["bytes"]]+as.numeric(object.size(value))
 }
}))
trace("dbExecute",where=asNamespace("DBI"),print=FALSE,exit=quote({
 if(isTRUE(metrics[["active"]])) metrics[["executes"]] = metrics[["executes"]]+1L
}))
records = list()
for (repetition in 0:3) {
 gc()
 metrics[["rows"]] = metrics[["bytes"]] = metrics[["reads"]] = metrics[["executes"]] = 0
 metrics[["active"]] = TRUE
 elapsed = system.time({fit = if(mode=="baseline") dbreg(y~x+z+f,conn=conn,table="dbreg_i_perf",strategy="moments",row_id="id") else
  dbreg::feols(y~x+z+i(f),conn=conn,table="dbreg_i_perf",row_id="id",lean=mode=="lean",notes=FALSE)})[["elapsed"]]
 metrics[["active"]] = FALSE
 stopifnot(max(abs(unname(coef(fit))-unname(coef(reference))))<1e-7,
  max(abs(unname(vcov(fit))-unname(vcov(reference))))<1e-9,nobs(fit)==1e6,
  abs(as.numeric(logLik(fit))-reference_ll)<1e-5,metrics[["rows"]]<=3L)
 if(mode=="lean") stopifnot(is.null(fit[["resource_state"]]),is.null(fit[["fitted_state"]][["prediction"]]))
 if(repetition>0) records[[repetition]] = data.frame(repetition=repetition,mode=mode,elapsed=elapsed,
  fit_read_calls=metrics[["reads"]],fit_execute_calls=metrics[["executes"]],
  max_rows=metrics[["rows"]],query_object_bytes=metrics[["bytes"]],model_object_bytes=as.numeric(object.size(fit)))
 close(fit)
 stopifnot(DBI::dbIsValid(conn),identical(sort(DBI::dbListTables(conn)),original))
 if(repetition>0) {
  status = readLines("/proc/self/status")
  rss = sub(".*:[[:space:]]*([0-9]+).*","\\1",status[grepl("^VmRSS:",status)])
  records[[repetition]][["rss_after_close_kib"]] = as.numeric(rss)
 }
 message(mode," repetition ",repetition,": ",elapsed,"s")
}
untrace("dbGetQuery",where=asNamespace("DBI"))
untrace("dbExecute",where=asNamespace("DBI"))
DBI::dbExecute(conn,"DROP VIEW dbreg_i_perf")
DBI::dbDisconnect(conn,shutdown=TRUE)
write.csv(do.call(rbind,records),paste0(root,".csv"),row.names=FALSE)
