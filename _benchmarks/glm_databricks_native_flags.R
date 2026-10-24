# Compare the native-query flags used by newer odbc releases, without installing
# anything or changing the user's connection configuration.
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("_benchmarks/teacher_cloglog_profile.R")
args=commandArgs(TRUE);out=args[[2]]
dir.create(file.path(out,"sql"),recursive=TRUE,showWarnings=FALSE)
log=new.env(parent=emptyenv());log[["out"]]=out;log[["id"]]=0L
log[["tag"]]=paste0("glm_native_flags_",Sys.getpid());log[["queries"]]=list()
saved=function(i) paste(readLines(file.path(args[[1]],"sql",paste0("teacher_integrated_174818_Databricks_integrated_1_",i,".sql"))),collapse="\n")
run=function() {
 conns=list()
 on.exit(for(conn in conns) DBI::dbDisconnect(conn@native))
 for(variant in c("control","native_flags")) {
  extra=if(variant=="native_flags") list(EnableNativeParameterizedQuery=1,PopulateParametersForNativeQuery=1) else list()
  native=do.call(DBI::dbConnect,c(list(drv=odbc::databricks(),
   driver=Sys.getenv("DATABRICKS_ODBC_DRIVER","/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
   httpPath=Sys.getenv("DATABRICKS_SQL_PATH"),useNativeQuery=TRUE,DefaultStringColumnLength="65536"),extra))
  log[["run"]]=paste0(variant,"_setup")
  conns[[variant]]=methods::new("teacher_profile_connection",native=native,log=log)
  DBI::dbExecute(native,"SET use_cached_result=false")
  DBI::dbExecute(native,paste("CREATE TEMP TABLE teacher_integrated_input_174818 AS",teacher_cloglog_sql(5e6)))
  DBI::dbExecute(conns[[variant]],saved(5))
 }
 results=list();reference=NULL
 for(variant in c("control","native_flags","native_flags","control")) {
  rep=sum(vapply(results,function(r) r[["variant"]]==variant,logical(1)))+1L
  log[["run"]]=paste0(variant,"_",rep)
  gc()
  timing=system.time({value=DBI::dbGetQuery(conns[[variant]],saved(10))})[["elapsed"]]
  order_rows=function(d) d[order(d[["__set"]],d[["__i"]],d[["__j"]],na.last=TRUE),]
  value=order_rows(value)
  if(is.null(reference)) reference=value
  stopifnot(nrow(reference)==nrow(value))
  for(name in setdiff(names(value),"__globals")) stopifnot(isTRUE(all.equal(as.numeric(value[[name]]),as.numeric(reference[[name]]),tolerance=1e-12)))
  stopifnot(isTRUE(all.equal(jsonlite::fromJSON(value[["__globals"]][value[["__set"]]==0]),
   jsonlite::fromJSON(reference[["__globals"]][reference[["__set"]]==0]),tolerance=1e-12)))
  row=data.frame(variant=variant,repetition=rep,elapsed=timing,rows=nrow(value),passed=TRUE)
  results[[length(results)+1L]]=row;print(row);flush.console()
  write.csv(do.call(rbind,results),file.path(out,"experiments.csv"),row.names=FALSE)
 }
}
run()
