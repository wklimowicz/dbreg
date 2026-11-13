library(DBI)
local({
  saved=options(dbreg.ap_solver="auto",dbreg.ap_method=NULL)
  on.exit(options(saved),add=TRUE)
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  set.seed(11117);n=600L
  d=data.frame(id=seq_len(n),x=rnorm(n),z=rnorm(n),g=sample(1:18,n,TRUE),h=sample(1:13,n,TRUE),k=sample(1:7,n,TRUE),w=runif(n,.2,2),o=rnorm(n))
  d[["endog"]]=.8*d[["z"]]+rnorm(n)
  d[["y"]]=d[["o"]]+.5*d[["x"]]+.6*d[["endog"]]+sin(d[["g"]])+cos(d[["h"]])+rnorm(n)
  DBI::dbWriteTable(conn,"fe_fixture",d)
  before=sort(DBI::dbListTables(conn))
  for(estimator in c("ols","iv"))for(three in c(FALSE,TRUE))for(kind in list("iid","hc1",~g)){
    fml=if(estimator=="ols"){if(three)y~x+z|g+h+k else y~x+z|g+h}else{if(three)y~x|g+h+k|endog~z else y~x|g+h|endog~z}
    fun=if(estimator=="ols")dbreg else dbivreg
    fit=fun(fml,conn=conn,table="fe_fixture",weights=~w,offset=~o,vcov=kind,strategy="demean",fe_max_iter=1L,fe_tol=1e-9)
    adjusted=d;adjusted[["y"]]=d[["y"]]-d[["o"]]
    ref=fixest::feols(fml,adjusted,weights=~w,vcov=if(identical(kind,"hc1"))"hetero"else kind,ssc=fixest::ssc(K.fixef="full"),fixef.rm="none",fixef.tol=1e-11,fixef.iter=10000L,notes=FALSE)
    ix=match(names(coef(fit)),sub("^fit_","",names(coef(ref))))
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)[ix]),tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)[ix,ix,drop=FALSE]),tolerance=1e-8)
    if(estimator=="ols"){
      X=model.matrix(if(three)~x+z+factor(g)+factor(h)+factor(k)else~x+z+factor(g)+factor(h),d)
      w=d[["w"]];y=d[["y"]]-d[["o"]]
      beta=as.numeric(qr.solve(X*sqrt(w),y*sqrt(w)))
      residual=y-as.numeric(X%*%beta);bread=solve(crossprod(X,w*X))
      if(identical(kind,"iid"))V=bread*sum(w*residual^2)/(n-ncol(X)) else {
        scores=X*(w*residual)
        if(identical(kind,"hc1")){meat=crossprod(scores);scale=n/(n-ncol(X))}else{
          meat=crossprod(rowsum(scores,d[["g"]]));G=length(unique(d[["g"]]));scale=G/(G-1)*(n-1)/(n-ncol(X))
        }
        V=scale*bread%*%meat%*%bread
      }
      positions=match(names(coef(fit)),colnames(X))
      expect_equal(as.numeric(coef(fit)),beta[positions],tolerance=1e-9)
      expect_equal(as.numeric(vcov(fit)),as.numeric(V[positions,positions,drop=FALSE]),tolerance=1e-9)
    }
    expect_identical(fit[["ap_method"]],"bounded_moments")
    expect_true(fit[["ap_iterations"]]<=1L)
    expect_true(fit[["ap_max_abs"]]<=1e-9);expect_true(fit[["ap_max_relative"]]<=1e-9)
    expect_true(fit[["ap_transfer_cells"]]<=getOption("dbreg.ap_moment_max_cells",16e6))
    expect_equal(sort(DBI::dbListTables(conn)),before)
    expect_identical(fit[["projection_diagnostics"]][[1L]][["status"]],"completed")
  }
  # Internal weight/residual spellings may be genuine source roles. Both
  # engines must preserve their typed identities and use the requested weights.
  DBI::dbExecute(conn,"CREATE TEMP TABLE collision_fe AS SELECT *,g AS __w,h AS x_tilde FROM fe_fixture")
  for(mode in c("moments","iterative")){
    prior=options(dbreg.ap_solver=mode)
    fit=dbreg(y~x+z|`__w`+x_tilde,conn=conn,table="collision_fe",weights=~w,offset=~o,vcov=~`__w`,strategy="demean",fe_tol=1e-9)
    options(prior)
    ref=fixest::feols(y~x+z|g+h,transform(d,y=y-o),weights=~w,vcov=~g,ssc=fixest::ssc(K.fixef="full"),fixef.rm="none",fixef.tol=1e-11,notes=FALSE)
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-8)
    expect_identical(sort(DBI::dbListTables(conn)),sort(c(before,"collision_fe")))
  }
  DBI::dbExecute(conn,"DROP TABLE collision_fe")
  # Dense ordinal dictionaries must preserve original BIGINT and separator keys.
  DBI::dbExecute(conn,"CREATE TEMP TABLE typed_fe AS SELECT * EXCLUDE(g,h), CAST('9007199254740992' AS BIGINT)+CAST(g AS BIGINT) AS g, CASE WHEN h%2=0 THEN 'a|b/' ELSE '' END || CAST(h AS VARCHAR) AS h FROM fe_fixture")
  fit=dbreg(y~x+z|g+h,conn=conn,table="typed_fe",weights=~w,offset=~o,strategy="demean",fe_max_iter=1L,fe_tol=1e-9)
  ref=fixest::feols(y~x+z|g+h,transform(d,y=y-o),weights=~w,fixef.rm="none",fixef.tol=1e-11,fixef.iter=10000L,notes=FALSE)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-8)
  expect_identical(fit[["ap_method"]],"bounded_moments")
  DBI::dbExecute(conn,"DROP TABLE typed_fe")
  # Disconnected components are structural null directions, not small pivots.
  dd=d;dd[["h"]]=dd[["g"]]+(dd[["h"]]%%3)*18
  DBI::dbWriteTable(conn,"disconnected_fe",dd)
  fit=dbreg(y~x|g+h,conn=conn,table="disconnected_fe",weights=~w,strategy="demean",fe_max_iter=1L,fe_tol=1e-9)
  ref=fixest::feols(y~x|g+h,dd,weights=~w,fixef.rm="none",fixef.tol=1e-11,fixef.iter=10000L,notes=FALSE)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ref)),tolerance=1e-8)
  expect_identical(fit[["ap_method"]],"bounded_moments")
  DBI::dbExecute(conn,"DROP TABLE disconnected_fe")
  # Uniform weight rescaling preserves coefficients and covariance. Raw weights
  # survive the projected relation so reported RSS retains its source units.
  ordinary=dbreg(y~x|g+h,conn=conn,table="fe_fixture",weights=~w,strategy="demean",vcov="hc1",fe_tol=1e-9)
  for(scale in c(1e-100,1e100)){
    DBI::dbExecute(conn,paste0("CREATE TEMP TABLE scaled_fe AS SELECT * EXCLUDE(w),w*",sprintf("%.17g",scale)," AS w FROM fe_fixture"))
    fit=dbreg(y~x|g+h,conn=conn,table="scaled_fe",weights=~w,strategy="demean",vcov="hc1",fe_tol=1e-9)
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(ordinary)),tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ordinary)),tolerance=1e-8)
    expect_equal(attr(vcov(fit),"rss")/scale,attr(vcov(ordinary),"rss"),tolerance=1e-8)
    DBI::dbExecute(conn,"DROP TABLE scaled_fe")
  }
  # Expected capacity rejections fall back without losing ownership records.
  caps=list(dbreg.duckdb_ap_max_entries=1,dbreg.duckdb_ap_max_edges=1,dbreg.ap_moment_max_cells=1,dbreg.projection_max_bytes=512,dbreg.projection_max_upload_bytes=1)
  for(name in names(caps)){
    prior=options(setNames(caps[name],name))
    fit=dbreg(y~x|g+h,conn=conn,table="fe_fixture",weights=~w,strategy="demean",fe_tol=1e-9)
    options(prior)
    expect_equal(as.numeric(coef(fit)),as.numeric(coef(ordinary)),tolerance=1e-8)
    expect_identical(fit[["ap_method"]],"row_tables")
    expect_identical(fit[["projection_diagnostics"]][[1L]][["status"]],"declined")
    expect_true(fit[["projection_diagnostics"]][[1L]][["reason"]] %in% c("system_limit","transfer_limit","upload_limit"))
    expect_identical(sort(DBI::dbListTables(conn)),before)
  }
  prior=options(dbreg.glm_fe_dense_max=0,dbreg.glm_fe_cross_max=1000)
  fit=dbreg(y~x|g+h,conn=conn,table="fe_fixture",weights=~w,strategy="demean",fe_tol=1e-9)
  options(prior)
  expect_identical(fit[["ap_solver"]],"schur_pcg")
  expect_identical(fit[["ap_preconditioner"]],"block_sgs")
  expect_true(fit[["ap_solver_iterations"]]>0L && fit[["ap_solver_iterations"]]<=2000L*8L)
  expect_equal(as.numeric(coef(fit)),as.numeric(coef(ordinary)),tolerance=1e-8)
  # Sparse Schur Cholesky is exercised through the real fitting consumer.
  # The exact DGP supplies an independent coefficient oracle on a sparse ring.
  ns=3200L
  ds=data.frame(g=rep(1:800,each=4),x=rnorm(ns))
  ds[["h"]]=(ds[["g"]]+rep(0:3,800))%%800+1
  ds[["y"]]=2*ds[["x"]]+sin(ds[["g"]])+cos(ds[["h"]])
  DBI::dbWriteTable(conn,"sparse_fe",ds)
  sparse=dbreg(y~x|g+h,conn=conn,table="sparse_fe",strategy="demean",fe_max_iter=4L,fe_tol=1e-8)
  expect_identical(sparse[["ap_solver"]],"sparse_cholesky")
  expect_equal(as.numeric(coef(sparse)),2,tolerance=1e-8)
  expect_true(sparse[["ap_max_relative"]]<=1e-8)
  DBI::dbExecute(conn,"DROP TABLE sparse_fe")
  expect_identical(sort(DBI::dbListTables(conn)),before)
  # Explicit iterative choices remain authoritative; a one-cycle budget fails.
  for(method in c("join","window")){
    prior=options(dbreg.ap_method=method)
    expect_error(dbreg(y~x|g+h,conn=conn,table="fe_fixture",weights=~w,strategy="demean",fe_max_iter=1L),"within 1 iterations")
    options(prior)
    expect_identical(sort(DBI::dbListTables(conn)),before)
  }
  prior=options(dbreg.ap_solver="iterative")
  expect_error(dbreg(y~x|g+h,conn=conn,table="fe_fixture",weights=~w,strategy="demean",fe_max_iter=1L),"within 1 iterations")
  options(prior)
  expect_identical(sort(DBI::dbListTables(conn)),before)
  # A rejected SQL row certificate cleans every private table. Unexpected
  # errors propagate through the policy, rather than becoming slow AP retries.
  inputs=list2env(list(conn=conn,from_statement="FROM fe_fixture",fe=c("g","h"),yvar="y",weights="w",cluster_var="g",verbose=FALSE,projection_limits=dbreg:::projection_limits()),parent=emptyenv())
  project=dbreg:::dbreg_moment_projections
  environment(project)=new.env(parent=environment(project))
  environment(project)[["fe_projection_certificate"]]=function(...)list(max_abs=Inf,max_relative=Inf)
  attempt=dbreg:::projection_attempt("bounded_moments",project(inputs,c("x","z"),c("x","z"),100L,1e-9))
  expect_identical(attempt[["reason"]],"certificate_failed")
  expect_identical(sort(DBI::dbListTables(conn)),before)
  environment(project)[["fe_projection_certificate"]]=function(...)stop("unexpected driver failure")
  expect_error(dbreg:::projection_attempt("bounded_moments",project(inputs,c("x","z"),c("x","z"),100L,1e-9)),"unexpected driver failure")
  expect_identical(sort(DBI::dbListTables(conn)),before)
  expect_true(DBI::dbIsValid(conn))
  for(bad in list(0,NA_real_,Inf,"small",c(1,2))){
    prior=options(dbreg.ap_moment_max_cells=bad)
    expect_error(dbreg(y~x|g+h,conn=conn,table="fe_fixture",strategy="demean"),"Projection limits")
    options(prior)
    expect_identical(sort(DBI::dbListTables(conn)),before)
  }
  # The row certificate uses weighted RMS, matching iterative AP, and
  # rejects missing values or a changed row count independently of the solve.
  cert=data.frame(f1=c(1,1,2,2),f2=c(1,2,1,2),w=1:4,v=c(-1,1,2,-2))
  DBI::dbWriteTable(conn,"certificate_fe",cert)
  measured=dbreg:::fe_projection_certificate(conn,"certificate_fe",c("f1","f2"),"v",4,2)
  means=unlist(lapply(cert[c("f1","f2")],function(group)tapply(seq_len(4),group,function(i)weighted.mean(cert[["v"]][i],cert[["w"]][i]))))
  expected=max(abs(means))
  expect_equal(measured[["max_abs"]],expected)
  expect_equal(measured[["max_relative"]],expected/sqrt(weighted.mean(cert[["v"]]^2,cert[["w"]])))
  expect_identical(dbreg:::fe_projection_certificate(conn,"certificate_fe",c("f1","f2"),"v",5,2)[["max_abs"]],Inf)
  DBI::dbExecute(conn,"UPDATE certificate_fe SET v=NULL WHERE w=1")
  expect_identical(dbreg:::fe_projection_certificate(conn,"certificate_fe",c("f1","f2"),"v",4,2)[["max_abs"]],Inf)
  DBI::dbExecute(conn,"DROP TABLE certificate_fe")
  prior=options(dbreg.projection_max_rows=100,dbreg.duckdb_ap_max_edges=1,
    dbreg.fe_dense_max=7,dbreg.glm_fe_dense_max=0)
  limits=dbreg:::projection_limits()
  expect_equal(limits[["summary_rows"]],100)
  expect_equal(limits[["dense_levels"]],7)
  options(prior)
  attr(conn,"connectionString")=list(options=list(UseNativeQuery=FALSE))
  expect_error(dbreg:::validate_native_temp_query(conn,list(is_databricks=TRUE)),"native query mode")
  attr(conn,"connectionString")=NULL
})
