# Exact same fixed working predictor, numerical checks before any fit changes.
# Rscript _benchmarks/glm_databricks_more_sql.R /tmp/teacher-cloglog-production /tmp/glm-more-sql
pkgload::load_all()
source("_benchmarks/teacher_cloglog_data.R")
source("_benchmarks/teacher_cloglog_profile.R")
source("_benchmarks/glm_databricks_more_helpers.R")
args=commandArgs(TRUE);out=args[[2]]
dir.create(file.path(out,"sql"),recursive=TRUE,showWarnings=FALSE)
log=new.env(parent=emptyenv());log[["out"]]=out;log[["id"]]=0L
log[["tag"]]=paste0("glm_more_",Sys.getpid());log[["run"]]="setup";log[["queries"]]=list()
saved=function(i) paste(readLines(file.path(args[[1]],"sql",paste0("teacher_integrated_174818_Databricks_integrated_1_",i,".sql"))),collapse="\n")
run=function() {
 native=DBI::dbConnect(odbc::databricks(),
  driver=Sys.getenv("DATABRICKS_ODBC_DRIVER","/opt/databricks/databricksodbc/lib/64/libdatabricksodbc64.so"),
  httpPath=Sys.getenv("DATABRICKS_SQL_PATH"),useNativeQuery=TRUE,DefaultStringColumnLength="65536")
 on.exit(DBI::dbDisconnect(native))
 conn=methods::new("teacher_profile_connection",native=native,log=log)
 DBI::dbExecute(native,"SET use_cached_result=false")
 DBI::dbExecute(native,paste("CREATE TEMP TABLE teacher_integrated_input_174818 AS",teacher_cloglog_sql(5e6)))
 DBI::dbExecute(conn,saved(5))
 working=glm_experiment_source(saved(10));fe=paste0("__dbglm_fe",1:6);x="dbreg_internal_x_0001__"
 values=dbreg:::dbglm_fe_value_terms(conn,c(x,"z"),"score");globals=glm_experiment_terms(conn,x)
 ordinary=function(source) dbreg:::dbglm_sql(conn,dbreg:::dbglm_databricks_moments_sql(source,fe,values,globals))
 original=ordinary(working)
 two=dbreg:::dbglm_sql(conn,glm_experiment_two_branches(working,fe,values,globals))
 fields=c("__i","__j","__h","__b1","__b2","__ex1")
 packed=paste0("SELECT __set,COUNT(*) AS __n,MAX(__globals) AS __globals,",
  "to_json(collect_list(array(",paste(fields,collapse=","),"))) AS __payload FROM (",original,
  ") p GROUP BY __set,COALESCE(__i,0)")
 results=list();reference=NULL
 for(rep in 1:2) for(variant in if(rep==1) c("control","two_branches","packed_fixed","materialized") else c("materialized","packed_fixed","two_branches","control")) {
  log[["run"]]=paste0(variant,"_",rep);gc();start=proc.time()[["elapsed"]]
  if(variant=="materialized") DBI::dbExecute(conn,paste("CREATE TEMP TABLE experiment_working AS SELECT",paste(c(fe,x,"eta","mu","w","z","score","dev"),collapse=","),"FROM",working))
  sql=switch(variant,control=original,two_branches=two,packed_fixed=packed,materialized=ordinary("experiment_working"))
  res=DBI::dbGetQuery(conn,sql)
  if(variant=="materialized") DBI::dbExecute(conn,"DROP TEMP TABLE experiment_working")
  query_elapsed=proc.time()[["elapsed"]]-start;decode=0;status="passed";numeric_error=0
  if(variant=="control" && is.null(reference)) reference=res
  decoded=res
  if(variant=="packed_fixed") {
   decode=system.time({error=tryCatch({
    chunks=lapply(seq_len(nrow(res)),function(i) {
      m=jsonlite::fromJSON(res[["__payload"]][[i]])
      stopifnot(nrow(m)==as.numeric(res[["__n"]][[i]]))
      d=as.data.frame(m);names(d)=fields;d[["__set"]]=res[["__set"]][[i]];d
    })
    decoded=do.call(rbind,chunks)
    decoded[["__globals"]]=NA_character_
    decoded[["__globals"]][decoded[["__set"]]==0]=res[["__globals"]][res[["__set"]]==0]
    NULL
   },error=function(e) conditionMessage(e))})[["elapsed"]]
   if(!is.null(error)) status=paste("rejected:",error)
  }
  if(status=="passed") {
   order_rows=function(d) d[order(d[["__set"]],d[["__i"]],d[["__j"]],na.last=TRUE),]
   a=order_rows(reference);b=order_rows(decoded)
   stopifnot(nrow(a)==nrow(b))
   for(name in c("__set",fields)) {
    av=as.numeric(a[[name]]);bv=as.numeric(b[[name]])
    stopifnot(identical(is.na(av),is.na(bv)),isTRUE(all.equal(av,bv,tolerance=1e-12)))
    numeric_error=max(numeric_error,max(abs(av-bv),na.rm=TRUE))
   }
   stopifnot(isTRUE(all.equal(jsonlite::fromJSON(a[["__globals"]][a[["__set"]]==0]),
     jsonlite::fromJSON(b[["__globals"]][b[["__set"]]==0]),tolerance=1e-12)))
  }
  row=data.frame(variant=variant,repetition=rep,query_seconds=query_elapsed,decode_seconds=decode,
   total_seconds=query_elapsed+decode,rows=nrow(res),result_bytes=as.numeric(object.size(res)),max_abs_error=numeric_error,status=status)
  results[[length(results)+1L]]=row;print(row);flush.console()
  write.csv(do.call(rbind,results),file.path(out,"experiments.csv"),row.names=FALSE)
 }
}
run()
