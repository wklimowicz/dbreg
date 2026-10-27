library(DBI)
library(dbreg)
source(system.file('tinytest/helper_glm_reference.R',package='dbreg'),local=TRUE)
source(system.file('tinytest/helper_glm_databricks.R',package='dbreg'),local=TRUE)
native=dbConnect(duckdb::duckdb(config=list(temp_directory='',threads='2')))
remote=dbglm_test_databricks_connection(native)
r = DBI::dbGetQuery(remote, paste0(
  "SELECT ", paste(dbreg:::sql_is_finite(remote, c("1.0", "'Infinity'::DOUBLE", "'NaN'::DOUBLE")),
    "AS", c("f", "i", "n"), collapse = ","),
  ",", dbreg:::sql_sum(remote, "x"), " AS s FROM (VALUES(1.0),(2.0)) t(x)"
))
expect_identical(as.logical(r[1,1:3]),c(TRUE,FALSE,FALSE))
expect_equal(r$s,3)
# Quoted identifiers and SQL string literals survive dialect translation.
r=DBI::dbGetQuery(remote,"SELECT \"FSUM\",'FSUM isfinite EXCLUDE VARCHAR' AS text FROM (SELECT 7 AS \"FSUM\") t")
expect_equal(r$FSUM,7)
expect_identical(r$text,'FSUM isfinite EXCLUDE VARCHAR')
v=data.frame(k=c("a'b",'FSUM'),x=c(1e-25,-1e25))
s=dbreg:::sql_values_frame(remote,v)
expect_equal(dbGetQuery(remote,s),v,check.attributes=FALSE)
# Extreme doubles and escaped keys survive JSON and the SQL parser exactly.
edge=data.frame(k=c("a'b",'x\\y','\n',intToUtf8(7),'é','東京'),
 x=c(.65313212622670436,.Machine$double.xmax,.Machine$double.xmin,1e-300,-1e300,.1))
expect_identical(dbGetQuery(remote,dbreg:::sql_values_frame(remote,edge))$k,edge$k)
expect_identical(dbGetQuery(remote,dbreg:::sql_values_frame(remote,edge))$x,edge$x)
expect_equal(dbGetQuery(remote,dbreg:::sql_values_frame(remote,v,format='values')),v,check.attributes=FALSE)
expect_null(dbreg:::sql_values_frame(remote,v,max_bytes=1))
set.seed(139)
d=data.frame(id=1:1200,g=rep(1:60,each=20),h=sample(1:10,1200,TRUE),x=rnorm(1200),z=rnorm(1200),a=sample(1:3,1200,TRUE),o=runif(1200,-.1,.1))
first=seq(1,1200,20);last=seq(20,1200,20)
for(family in list(binomial('logit'),binomial('probit'),binomial('cloglog'),poisson())) {
 eta=-.3+.2*d$x-.15*d$z+.1*sin(d$g)+d$o
 d$y=if(family$family=='poisson')rpois(1200,exp(eta)) else rbinom(1200,1,family$linkinv(eta))
 if(family$family=='poisson')d$y[first]=1 else {d$y[first]=0;d$y[last]=1}
 dbWriteTable(native,'input',d,overwrite=TRUE)
 for(kind in c('iid','hc1','cluster')) {
  vc=if(kind=='cluster')~h else kind
  ref=dbglm_reference_feglm(y~x+z|g+h,data=d,family=family,weights=~a,offset=~o,vcov=vc)
  m=dbglm(y~x+z|g+h,conn=remote,table='input',family=family,weights='a',offset='o',vcov=vc,row_id='id')
  expect_true(max(abs(coef(m)-coef(ref)))<1e-6)
  expect_true(max(abs(vcov(m)-vcov(ref)))<1e-6)
  expect_true(abs(as.numeric(logLik(m)-logLik(ref)))/nrow(d)<1e-9)
  expect_identical(m$projection_method,'moments')
  expect_identical(dbListTables(native),'input')
  # FE-level PCG remains the two-way path when the dense solve is capped.
  saved=options(dbreg.glm_fe_dense_max=0)
  m=dbglm(y~x+z|g+h,conn=remote,table='input',family=family,weights='a',offset='o',vcov=vc,row_id='id')
  options(saved)
  expect_true(max(abs(coef(m)-coef(ref)))<1e-6)
  expect_true(max(abs(vcov(m)-vcov(ref)))<1e-6)
  expect_identical(m$projection_method,'bounded_pcg')
  expect_identical(dbListTables(native),'input')
 }
 # The original Spark window path now has the same GLM certificate.
 for(chunk in c(1L,3L)) {
  saved=options(dbreg.spark_ap_chunk_cycles=chunk)
  m=dbglm(y~x+z|g+h,conn=remote,table='input',family=family,control=dbglm_control(projection_method='sql'))
  options(saved)
  ref=dbglm_reference_feglm(y~x+z|g+h,data=d,family=family)
  expect_true(max(abs(coef(m)-coef(ref)))<1e-6)
  expect_true(max(abs(vcov(m)-vcov(ref)))<1e-6)
  expect_true(m$projection_residual<=m$control$fe_tol)
  expect_identical(m$projection_method,'window')
  expect_identical(dbListTables(native),'input')
 }
}
# Cap fallback must release the FE-summary views as well as physical tables.
saved=options(dbreg.duckdb_ap_max_edges=1)
fallback=dbglm(y~x+z|g+h,conn=remote,table='input',family=poisson())
options(saved)
ref=dbglm_reference_feglm(y~x+z|g+h,data=d,family=poisson())
expect_identical(fallback$projection_method,'window')
expect_true(max(abs(coef(fallback)-coef(ref)))<1e-6)
expect_identical(dbListTables(native),'input')
# Multi-way (>2 FE) models on Databricks use the single-pass moments projector.
d$k1=sample(1:6,1200,TRUE);d$k2=sample(1:4,1200,TRUE)
d$y=as.integer(d$y>0);d$y[first]=0;d$y[last]=1
dbWriteTable(native,'input',d,overwrite=TRUE)
for(kind in c('iid','hc1','cluster')) {
 vc=if(kind=='cluster')~h else kind
 ref_m=dbglm_reference_feglm(y~x+z|g+h+k1+k2,data=d,family=binomial('cloglog'),weights=~a,offset=~o,vcov=vc)
 m_m=dbglm(y~x+z|g+h+k1+k2,conn=remote,table='input',family=binomial('cloglog'),weights='a',offset='o',vcov=vc,row_id='id')
 expect_identical(m_m$projection_method,'moments')
 expect_true(max(abs(coef(m_m)-coef(ref_m)))<1e-6)
 expect_true(max(abs(vcov(m_m)-vcov(ref_m)))<1e-6)
 expect_identical(dbListTables(native),'input')
}
saved=options(dbreg.duckdb_ap_max_edges=1)
fallback_m=dbglm(y~x+z|g+h+k1+k2,conn=remote,table='input',family=binomial('cloglog'))
options(saved)
expect_identical(fallback_m$projection_method,'window')
expect_identical(dbListTables(native),'input')
# A failed fit removes its own state and preserves the source.
expect_error(dbglm(y~x+z|g+h,conn=remote,table='input',family=poisson(),control=dbglm_control(glm_max_iter=1)),'converge')
expect_identical(dbListTables(native),'input')
dbDisconnect(native,shutdown=TRUE)
