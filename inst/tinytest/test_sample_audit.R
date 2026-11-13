library(DBI)
local({
  set.seed(11104);n=160L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=as.numeric(rep(1:8,20)),w=runif(n,.2,2),cl=as.numeric(rep(1:10,16)),o=rnorm(n),keep=rep(c(TRUE,TRUE,FALSE,TRUE),40),u=rnorm(n),f=rep(c("NaN","a","b","c"),40))
  d[["endog"]]=.7*d[["z"]]+rnorm(n);d[["y"]]=d[["o"]]+.5*d[["x"]]+.7*d[["endog"]]+rnorm(n)
  for(name in c("x","y","g","w","cl","endog","z"))d[[name]][match(name,c("x","y","g","w","cl","endog","z"))]=NA
  d[["w"]][8]=0;d[["o"]][9:11]=c(NA,Inf,NaN);d[["keep"]][12]=NA;d[["u"]][13:14]=c(NA,NaN)
  subset=d[["keep"]] & d[["u"]]>0;subset[is.na(subset)]=FALSE
  for(estimator in c("ols","iv"))for(strategy in if(estimator=="ols")c("moments","compress","demean","mundlak")else c("moments","compress","demean")){
    has_fe=strategy %in% c("demean","mundlak")
    formula=if(estimator=="ols"){if(has_fe)y~x+factor(f)|g else y~x+factor(f)}else{if(has_fe)y~x|g|endog~z else y~x|endog~z}
    fun=if(estimator=="ols")dbreg else dbivreg
    args=list(fml=formula,data=d,weights="w",offset="o",subset=~keep & u>0,strategy=strategy,vcov=~cl)
    fit=do.call(fun,c(args,list(diagnostics=TRUE)));plain=do.call(fun,args)
    audit=fit[["sample_diagnostics"]]
    expected_roles=if(estimator=="ols")list(outcome="y",regressors=c("x","f"))else list(outcome="y",regressors="x",endogenous="endog",instruments="z")
    expected_roles=c(expected_roles,if(has_fe)list(fe="g"),list(weights="w",cluster="cl",offset="o"))
    fields=unlist(expected_roles,use.names=FALSE)
    selected=subset & complete.cases(d[unique(fields)]) & is.finite(d[["y"]]-d[["o"]]) & d[["w"]]>0
    expect_identical(audit[["n_input"]],as.numeric(n))
    expect_identical(fit[["nobs_input"]],as.numeric(n))
    expect_equal(audit[["n_selected"]],sum(subset))
    expect_equal(audit[["n_fitted"]],sum(selected));expect_equal(nobs(fit),sum(selected))
    expect_equal(audit[["n_dropped"]],n-sum(selected))
    expect_equal(audit[["n_nonpositive_weight"]],sum(d[["w"]]<=0,na.rm=TRUE))
    expect_equal(audit[["n_invalid_offset"]],sum(!is.finite(d[["y"]]) | !is.finite(d[["o"]]) | !is.finite(d[["y"]]-d[["o"]])))
    expect_identical(audit[["detail"]][["variable"]],fields)
    expect_identical(audit[["detail"]][["role"]],rep(names(expected_roles),lengths(expected_roles)))
    expect_equal(audit[["detail"]][["n_missing"]],vapply(fields,function(name)sum(is.na(d[[name]])),numeric(1)),check.attributes=FALSE)
    expect_equal(coef(fit,fe=TRUE),coef(plain,fe=TRUE),tolerance=1e-12)
    expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(vcov(plain,fe=TRUE)),tolerance=1e-12)
    expect_true(is.na(plain[["nobs_input"]]));expect_true(is.null(plain[["sample_diagnostics"]]))
    if(estimator=="iv")expect_equal(fit[["diagnostics"]][["first_stage_f"]],plain[["diagnostics"]][["first_stage_f"]],tolerance=1e-12)
  }
  # Source-role maps refer to the original, public names rather than the target.
  names(d)[match(c("y","o","w","keep"),names(d))]=c("outcome value","offset value","case weight","keep row")
  fit=dbreg(`outcome value`~x+I(x^2),data=d,offset=~`offset value`,weights="case weight",subset=~`keep row` & u>0,diagnostics=TRUE,strategy="moments")
  expect_true(all(c("outcome value","offset value","case weight") %in% fit[["sample_diagnostics"]][["detail"]][["variable"]]))
  expect_identical(fit[["yvar"]],"outcome value")
  valid=subset & complete.cases(d[c("outcome value","x","offset value","case weight")]) & is.finite(d[["outcome value"]]-d[["offset value"]]) & d[["case weight"]]>0
  reference=lm(I(`outcome value`-`offset value`)~x+I(x^2),d[valid,],weights=`case weight`)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(reference)),tolerance=1e-8)
  expect_equal(predict(fit,d[valid,]),as.numeric(predict(reference,d[valid,]))+d[["offset value"]][valid],tolerance=1e-8)
  names(d)[match(c("g","endog","z"),names(d))]=c("fe value","endog value","instrument value")
  iv=dbivreg(`outcome value`~x|`fe value`|`endog value`~`instrument value`,data=d,
    offset="offset value",weights="case weight",subset=~`keep row` & u>0,diagnostics=TRUE,strategy="demean")
  valid=subset & complete.cases(d[c("outcome value","x","fe value","endog value","instrument value","offset value","case weight")]) & is.finite(d[["outcome value"]]-d[["offset value"]]) & d[["case weight"]]>0
  rows=d[valid,];w=rows[["case weight"]];F=model.matrix(~factor(`fe value`),rows)
  U=cbind(as.matrix(rows[c("x","endog value","instrument value")]),rows[["outcome value"]]-rows[["offset value"]])
  projected=U-qr.fitted(qr(F*sqrt(w)),U*sqrt(w))/sqrt(w)
  X=projected[,1:2];Z=projected[,c(1,3)];y=projected[,4]
  A=crossprod(X,w*Z)%*%solve(crossprod(Z,w*Z));beta=as.numeric(solve(A%*%crossprod(Z,w*X),A%*%crossprod(Z,w*y)))
  expect_equal(as.numeric(coef(iv)),beta,tolerance=1e-8)
  expect_identical(iv[["yvar"]],"outcome value")
  expect_true(all(c("outcome value","fe value","endog value","instrument value") %in% iv[["sample_diagnostics"]][["detail"]][["variable"]]))
  for(bad in list(NA,1,"yes",c(TRUE,FALSE),logical())){
    expect_error(dbreg(`outcome value`~x,data=d,diagnostics=bad),"TRUE or FALSE")
    expect_error(dbivreg(`outcome value`~x|`endog value`~`instrument value`,data=d,diagnostics=bad),"TRUE or FALSE")
  }
})

local({
  conn=DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11105);n=80L
  d=data.frame(id=seq_len(n),x=rnorm(n),g=as.numeric(rep(1:8,10)),w=runif(n,.2,2),cl=as.numeric(rep(1:10,8)),o=rnorm(n),keep=rep(c(TRUE,FALSE),40),f=rep(c("NaN","a"),40))
  d[["y"]]=d[["o"]]+.5*d[["x"]]+rnorm(n)
  DBI::dbWriteTable(conn,"audit_source",d)
  for(name in c("x","y","g","w","cl","o"))DBI::dbExecute(conn,paste0("UPDATE audit_source SET ",name,"=CAST('NaN' AS DOUBLE) WHERE id=",match(name,c("x","y","g","w","cl","o"))))
  fit=dbreg(y~x+factor(f)|g,conn=conn,table="audit_source",weights="w",offset="o",vcov=~cl,strategy="demean",diagnostics=TRUE)
  audit=fit[["sample_diagnostics"]]
  expect_equal(audit[["detail"]][["n_missing"]],c(1,1,0,1,1,1,1))
  expect_equal(audit[["n_fitted"]],n-6);expect_equal(audit[["n_invalid_offset"]],2)
  expect_equal(audit[["n_selected"]],n)
  replay=DBI::dbGetQuery(conn,audit[["query_string"]])
  expect_equal(as.numeric(replay[[audit[["count_aliases"]][["n_input"]]]]),n)
  # Distinct source fields and generated count aliases share similar names.
  names(d)[names(d)=="x"]="dbreg_audit_n_input"
  fit=dbreg(y~dbreg_audit_n_input,data=d,subset=~dbreg_audit_n_input>0,diagnostics=TRUE,strategy="moments")
  expect_equal(fit[["sample_diagnostics"]][["n_selected"]],sum(d[["dbreg_audit_n_input"]]>0))
  expect_identical(fit[["sample_diagnostics"]][["detail"]][["variable"]],c("y","dbreg_audit_n_input"))
  # Counts remain zero on an empty typed source, even without fitted results.
  DBI::dbWriteTable(conn,"empty_audit",d[FALSE,])
  inputs=new.env(parent=emptyenv());inputs[["diagnostics"]]=TRUE;inputs[["conn"]]=conn
  source=dbreg:::safe_source(conn,"FROM empty_audit",c("y","dbreg_audit_n_input"))
  audit=dbreg:::sample_diagnostics(inputs,source[["from_statement"]],list(outcome="y",regressors=source[["to_safe"]]("dbreg_audit_n_input")),source)
  expect_equal(c(audit[["n_input"]],audit[["n_selected"]],audit[["n_nonpositive_weight"]],audit[["n_invalid_offset"]]),rep(0,4))
  expect_equal(audit[["detail"]][["n_missing"]],c(0,0))
})

local({
  # This DBI wrapper observes query work while DuckDB executes real SQL.
  connection_class="dbreg_audit_query_connection"
  methods::setClass(connection_class,contains="DBIConnection",slots=c(native="ANY",state="environment"))
  methods::setMethod("dbGetInfo",connection_class,function(dbObj,...)list(dbms.name="DuckDB"))
  methods::setMethod("dbQuoteIdentifier",c(connection_class,"character"),function(conn,x,...)DBI::dbQuoteIdentifier(conn@native,x,...))
  methods::setMethod("dbGetQuery",c(connection_class,"character"),function(conn,statement,...){
    conn@state[["queries"]]=c(conn@state[["queries"]],statement)
    DBI::dbGetQuery(conn@native,statement,...)
  })
  methods::setMethod("dbSendQuery",c(connection_class,"character"),function(conn,statement,...){
    conn@state[["queries"]]=c(conn@state[["queries"]],statement)
    DBI::dbSendQuery(conn@native,statement,...)
  })
  native=DBI::dbConnect(duckdb::duckdb())
  on.exit({
    DBI::dbDisconnect(native,shutdown=TRUE)
    methods::removeMethod("dbGetInfo",connection_class)
    methods::removeMethod("dbQuoteIdentifier",c(connection_class,"character"))
    methods::removeMethod("dbGetQuery",c(connection_class,"character"))
    methods::removeMethod("dbSendQuery",c(connection_class,"character"))
    methods::removeClass(connection_class)
  })
  set.seed(11106);d=data.frame(x=rnorm(100),y=rnorm(100),z=rnorm(100),endog=rnorm(100))
  DBI::dbWriteTable(native,"audit_queries",d)
  state=new.env(parent=emptyenv());state[["queries"]]=character()
  conn=methods::new(connection_class,native=native,state=state)
  for(estimator in c("ols","iv")){
    fun=if(estimator=="ols")dbreg else dbivreg
    formula=if(estimator=="ols")y~x else y~x|endog~z
    state[["queries"]]=character();plain=fun(formula,conn=conn,table="audit_queries",strategy="moments")
    baseline=length(state[["queries"]])
    state[["queries"]]=character();audited=fun(formula,conn=conn,table="audit_queries",strategy="moments",diagnostics=TRUE)
    expect_equal(length(state[["queries"]]),baseline+2)
    expect_equal(coef(audited),coef(plain),tolerance=1e-12)
    for(mode in c("sql_only","data_only")){
      args=list(fml=formula,conn=conn,table="audit_queries",strategy="moments");args[[mode]]=TRUE
      state[["queries"]]=character();invisible(capture.output({first=do.call(fun,args)}))
      baseline=length(state[["queries"]])
      state[["queries"]]=character();invisible(capture.output({second=do.call(fun,c(args,list(diagnostics=TRUE)))}))
      expect_equal(length(state[["queries"]]),baseline)
      expect_equal(second,first,tolerance=1e-12)
    }
  }
})
