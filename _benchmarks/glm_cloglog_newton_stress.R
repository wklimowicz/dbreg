# Numerical screen for the experimental solver; no production defaults change.
pkgload::load_all()
source("_benchmarks/glm_databricks_more_helpers.R")
source("inst/tinytest/helper_glm_reference.R")
out=commandArgs(TRUE)[[1]];dir.create(out,recursive=TRUE,showWarnings=FALSE)
conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2",temp_directory="")))
set.seed(8942)
n=4500L
d=data.frame(f1=rep(1:18,length.out=n),f2=sample(1:10,n,TRUE),f3=sample(1:7,n,TRUE),
 f4=sample(1:5,n,TRUE),f5=rbinom(n,1,.7),f6=rbinom(n,1,.3),u=rnorm(n),z=rnorm(n),t=runif(n,-1,1))
d[["f2"]][1:3000]=(d[["f1"]][1:3000]-1L)%%10L+1L
cases=list(
 rare=list(rate=.02,weights=FALSE,offset=FALSE,shift=0,fe=6),
 very_rare_weighted=list(rate=.005,weights=TRUE,offset=TRUE,shift=0,fe=6),
 common_weighted=list(rate=.35,weights=TRUE,offset=TRUE,shift=0,fe=6),
 two_fe=list(rate=.6,weights=FALSE,offset=TRUE,shift=0,fe=2),
 shifted=list(rate=.1,weights=TRUE,offset=TRUE,shift=1e6,fe=6),
 varying_slope=list(rate=.15,weights=TRUE,offset=TRUE,shift=0,fe=3))
results=list()
for(name in names(cases)) {
 cfg=cases[[name]];data=d
 data[["x"]]=10.8+cfg[["shift"]]+.5*sin(data[["f1"]])+.3*data[["u"]]
 data[["a"]]=if(cfg[["weights"]]) sample(1:3,n,TRUE) else 1
 data[["o"]]=if(cfg[["offset"]]) runif(n,-.2,.2) else 0
 eta=log(-log(1-cfg[["rate"]]))+.25*data[["u"]]-.15*data[["z"]]+.25*sin(data[["f1"]])+data[["o"]]
 if(name=="varying_slope") eta=eta+.2*cos(data[["f1"]])*data[["t"]]
 data[["y"]]=rbinom(n,1,-expm1(-exp(eta)))
 formula=switch(as.character(cfg[["fe"]]),
  "2"=y~x+z|f1+f2,"3"=y~x+z|f1[t]+f2+f3,
  y~x+z|f1+f2+f3+f4+f5+f6)
 reference=dbglm_reference_feglm(formula,data,family=binomial("cloglog"),weights=~a,offset=~o,vcov=~f2,fixef.rm="perfect_fit")
 for(variant in c("control","adaptive_newton")) {
  experiment=if(variant=="adaptive_newton") glm_experiment_newton("adaptive") else NULL
  timing=system.time({answer=tryCatch(list(fit=dbglm(formula,data=data,conn=conn,family=binomial("cloglog"),
   weights="a",offset="o",vcov=~f2)),error=function(e) list(error=conditionMessage(e)),
   finally=if(!is.null(experiment)) experiment[["restore"]]())})
  fit=answer[["fit"]]
  passed=!is.null(fit)
  if(passed) passed=identical(names(coef(fit)),names(coef(reference))) &&
    all(abs(coef(fit)-coef(reference))<=1e-7+1e-6*abs(coef(reference))) &&
    all(abs(vcov(fit)-vcov(reference))<=1e-8+1e-5*abs(vcov(reference))) &&
    nobs(fit)==nobs(reference) && fit[["score_max"]]<=fit[["control"]][["score_tol"]] &&
    all(diff(fit[["deviance_history"]])<=64*.Machine$double.eps*pmax(1,head(fit[["deviance_history"]],-1))) &&
    (is.null(experiment)||experiment[["state"]][["last_fisher"]])
  row=data.frame(case=name,variant=variant,elapsed=timing[["elapsed"]],
   iterations=if(is.null(fit)) NA else fit[["iterations"]],
   halvings=if(is.null(fit)) NA else fit[["step_halvings"]],
   coef_error=if(is.null(fit)) NA else max(abs(coef(fit)-coef(reference))),
   vcov_error=if(is.null(fit)) NA else max(abs(vcov(fit)-vcov(reference))),
   nobs=if(is.null(fit)) NA else nobs(fit),passed=passed,
   error=if(is.null(answer[["error"]])) "" else answer[["error"]])
  results[[length(results)+1L]]=row;print(row);flush.console()
  write.csv(do.call(rbind,results),file.path(out,"cases.csv"),row.names=FALSE)
 }
}
DBI::dbDisconnect(conn,shutdown=TRUE)
