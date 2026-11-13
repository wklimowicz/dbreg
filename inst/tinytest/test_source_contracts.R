local({
  conn=DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11108);n=120L
  d=data.frame(x=rnorm(n),z=rnorm(n),w=runif(n,.2,2),g=rep(1:10,12),cl=rep(1:12,10),o=rnorm(n),keep=rep(c(TRUE,FALSE),60))
  d[["endog"]]=.7*d[["z"]]+rnorm(n);d[["y"]]=d[["o"]]+.5*d[["x"]]+.7*d[["endog"]]+rnorm(n)
  DBI::dbWriteTable(conn,"chosen_source",d)
  ignored=d[1:3,];ignored[["w"]]=-1
  for(estimator in c("ols","iv"))for(strategy in if(estimator=="ols")c("moments","compress","demean","mundlak")else c("moments","compress","demean")){
    has_fe=strategy %in% c("demean","mundlak")
    formula=if(estimator=="ols"){if(has_fe)y~x+z|g else y~x+z}else{if(has_fe)y~x|g|endog~z else y~x|endog~z}
    fun=if(estimator=="ols")dbreg else dbivreg
    args=list(fml=formula,conn=conn,table="chosen_source",weights=~w,offset=~o,subset=~keep,strategy=strategy,vcov=~cl,diagnostics=TRUE)
    fit=do.call(fun,c(args,list(data=ignored)))
    plain=do.call(fun,args)
    expect_equal(coef(fit,fe=TRUE),coef(plain,fe=TRUE),tolerance=1e-12)
    expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(vcov(plain,fe=TRUE)),tolerance=1e-12)
    expect_equal(nobs(fit),sum(d[["keep"]]))
    expect_equal(fit[["sample_diagnostics"]][["n_input"]],n)
    expect_identical(fit[["weights"]],"w")
  }
  for(estimator in c("ols","iv")){
    fun=if(estimator=="ols")dbreg else dbivreg
    formula=if(estimator=="ols")y~x else y~x|endog~z
    fit=fun(formula,conn=conn,table="chosen_source",data=stop("ignored data forced"),path=stop("ignored path forced"),weights=~w,strategy="moments",diagnostics=TRUE)
    expect_equal(nobs(fit),n);expect_equal(fit[["sample_diagnostics"]][["n_input"]],n)
    from_data=fun(formula,data=d,path=stop("ignored path forced"),weights=~w,strategy="moments")
    expect_equal(coef(from_data),coef(fit),tolerance=1e-9)
    if(estimator=="ols"){
      reference=lm(y~x,d,weights=w)
      expect_equal(as.numeric(coef(fit)),as.numeric(coef(reference)),tolerance=1e-9)
    }else{
      X=model.matrix(~x+endog,d);Z=model.matrix(~x+z,d);w=d[["w"]]
      P=crossprod(X,w*Z)%*%solve(crossprod(Z,w*Z))
      beta=as.numeric(solve(P%*%crossprod(Z,w*X),P%*%crossprod(Z,w*d[["y"]])))
      expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-9)
    }
    DBI::dbExecute(conn,"UPDATE chosen_source SET w=-1 WHERE rowid=0")
    expect_error(fun(formula,conn=conn,table="chosen_source",data=d,weights=~w,strategy="moments"),"non-negative")
    DBI::dbExecute(conn,"UPDATE chosen_source SET w=CAST('Infinity' AS DOUBLE) WHERE rowid=0")
    expect_error(fun(formula,conn=conn,table="chosen_source",data=ignored,weights="w",strategy="moments"),"finite")
    DBI::dbWriteTable(conn,"chosen_source",d,overwrite=TRUE)
    automatic=fun(formula,conn=conn,table="chosen_source",data=ignored,weights=~w,strategy="auto")
    reference=fun(formula,conn=conn,table="chosen_source",weights="w",strategy="auto")
    expect_identical(automatic[["strategy"]],reference[["strategy"]])
    expect_equal(coef(automatic),coef(reference),tolerance=1e-12)
    for(weight in list(~w+x,~log(w),~1,y~w)){
      expect_error(fun(formula,conn=conn,table="chosen_source",weights=weight),"formulas must name one column")
    }
    expect_true(DBI::dbIsValid(conn))
    expect_equal(DBI::dbListTables(conn),"chosen_source")
  }
  # The inherited source contract also governs binscatter inputs.
  bins=dbbinsreg(y~x,conn=conn,table="chosen_source",data=stop("ignored data forced"),path=stop("ignored path forced"),nbins=6,ci=FALSE,plot=FALSE)
  reference_bins=dbbinsreg(y~x,conn=conn,table="chosen_source",nbins=6,ci=FALSE,plot=FALSE)
  expect_equal(bins[["points"]],reference_bins[["points"]],tolerance=1e-12)
  # Resolve quoted column formulas without evaluating their environments.
  names(d)[names(d)=="w"]="case weight"
  fit=dbreg(y~x,data=d,weights=~`case weight`,strategy="moments")
  reference=lm(y~x,d,weights=`case weight`)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(reference)),tolerance=1e-9)
  expect_identical(fit[["weights"]],"case weight")
  expect_error(dbreg(y~x,data=d,weights=~missing),"not found")
  # A selected file still works after lower-priority promises are normalized.
  path=tempfile(fileext=".csv");on.exit(unlink(path),add=TRUE)
  write.csv(d,path,row.names=FALSE)
  fit=dbreg(y~x,path=path,weights=~`case weight`,strategy="moments")
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(reference)),tolerance=1e-8)
})

if(requireNamespace("dplyr",quietly=TRUE)&&requireNamespace("dbplyr",quietly=TRUE))local({
  conn=DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11109);n=100L
  d=data.frame(x=rnorm(n),z=rnorm(n),w=runif(n,.2,2))
  d[["endog"]]=.7*d[["z"]]+rnorm(n);d[["y"]]=.5*d[["x"]]+.7*d[["endog"]]+rnorm(n)
  DBI::dbWriteTable(conn,"lazy_source",d)
  lazy=dplyr::filter(dplyr::tbl(conn,"lazy_source"),x>0)
  for(estimator in c("ols","iv")){
    fun=if(estimator=="ols")dbreg else dbivreg
    formula=if(estimator=="ols")y~x else y~x|endog~z
    fit=fun(formula,data=lazy,path=stop("ignored path forced"),weights=~w,diagnostics=TRUE,strategy="moments")
    reference=fun(formula,table=lazy,weights="w",diagnostics=TRUE,strategy="moments")
    expect_equal(coef(fit),coef(reference),tolerance=1e-12)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(reference)),tolerance=1e-12)
    expect_equal(fit[["nobs_input"]],sum(d[["x"]]>0))
    expect_true(DBI::dbIsValid(conn))
  }
  bins=dbbinsreg(y~x,data=lazy,path=stop("ignored path forced"),nbins=5,ci=FALSE,plot=FALSE)
  reference_bins=dbbinsreg(y~x,table=lazy,nbins=5,ci=FALSE,plot=FALSE)
  expect_equal(bins[["points"]],reference_bins[["points"]],tolerance=1e-12)
  expect_true(DBI::dbIsValid(conn))

})
