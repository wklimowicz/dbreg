pkgload::load_all()
source("_benchmarks/glm_databricks_more_helpers.R")
input=commandArgs(TRUE)[[1]];out=commandArgs(TRUE)[[2]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="4",memory_limit="4GB",temp_directory=file.path(out,"tmp"))))
DBI::dbExecute(conn,paste0("CREATE TEMP TABLE teachers AS SELECT * FROM read_parquet(",DBI::dbQuoteString(conn,file.path(input,"teachers.parquet")),")"))
reference=readRDS(file.path(input,"fixest_strict_1.rds"))
results=list()
for(switch_at in c(1L,5L,6L,7L)) {
  cat("START switch_at",switch_at,"\n");flush.console()
  experiment=glm_experiment_newton(switch_at)
  timing=system.time({answer=tryCatch(list(fit=dbglm(wastage~log_outside_salary|spine_band_year+district+Experience+AGE+SEX+FTPTWK,
    conn=conn,table="teachers",family=binomial("cloglog"),vcov=~district,verbose=TRUE)),
    error=function(e) list(error=conditionMessage(e)),finally=experiment[["restore"]]())})
  fit=answer[["fit"]];error=answer[["error"]]
  row=data.frame(switch_at=switch_at,elapsed=timing[["elapsed"]],status=if(is.null(error)) "passed" else error,
    iterations=if(is.null(fit)) NA else fit[["iterations"]],
    coef_error=if(is.null(fit)) NA else max(abs(coef(fit)-reference[["coefficients"]])),
    se_error=if(is.null(fit)) NA else max(abs(sqrt(diag(vcov(fit)))-sqrt(diag(reference[["vcov"]])))),
    score=if(is.null(fit)) NA else fit[["score_max"]],last_fisher=experiment[["state"]][["last_fisher"]])
  print(row);flush.console();results[[length(results)+1L]]=row
  write.csv(do.call(rbind,results),file.path(out,"fits.csv"),row.names=FALSE)
  if(!is.null(fit)) stopifnot(row[["coef_error"]]<1e-8,row[["se_error"]]<1e-7,row[["last_fisher"]])
}
DBI::dbDisconnect(conn,shutdown=TRUE)
