library(DBI)
library(dbreg)
source(system.file('tinytest/helper_glm_reference.R',package='dbreg'),local=TRUE)
# Execute generated Databricks SQL on a RAM-only SQL fixture. Only the snapshot
# ID, star-exclusion and DROP syntax are adapted; this is a portability test, not remote timing.
methods::setClass('dbglm_sql_fixture',contains='DBIConnection',slots=c(native='ANY'))
methods::setMethod('dbGetInfo','dbglm_sql_fixture',function(dbObj,...)list(dbms.name='Databricks SQL'))
methods::setMethod('dbIsValid','dbglm_sql_fixture',function(dbObj,...)dbIsValid(dbObj@native))
methods::setMethod('dbQuoteIdentifier',c('dbglm_sql_fixture','character'),function(conn,x,...)dbQuoteIdentifier(conn@native,x,...))
methods::setMethod('dbQuoteString',c('dbglm_sql_fixture','character'),function(conn,x,...)dbQuoteString(conn@native,x,...))
methods::setMethod('dbExistsTable',c('dbglm_sql_fixture','character'),function(conn,name,...)dbExistsTable(conn@native,name,...))
fixture_sql=function(statement) {
 # Spark's JSON reader is emulated with R in this SQL fixture; actual Spark
 # decoding, extreme doubles and key escaping are checked remotely as well.
 while(grepl(" FROM explode(from_json(",statement,fixed=TRUE)) {
  m_pos=regexpr("SELECT CAST(item[",statement,fixed=TRUE)[1]
  sub_stmt=substr(statement,m_pos,nchar(statement))
  end_pos=regexpr("AS lookup(item)",sub_stmt,fixed=TRUE)[1]+nchar("AS lookup(item)")-1L
  block=substr(sub_stmt,1,end_pos)
  encoded=regmatches(block,regexec("unbase64\\('([^']*)'\\)",block,perl=TRUE))[[1]][2]
  rows=jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(encoded)),simplifyVector=FALSE)
  fields=regmatches(block,gregexpr('CAST\\(item\\[[0-9]+\\] AS (?:DOUBLE|STRING)\\) AS (?:"[^\"]+"|[A-Za-z_][A-Za-z0-9_]*)',block,perl=TRUE))[[1]]
  frame=as.data.frame(lapply(seq_along(fields),function(j) {
   if(grepl(' AS DOUBLE)',fields[j],fixed=TRUE)) vapply(rows,function(r)as.numeric(r[[j]]),numeric(1)) else
     vapply(rows,function(r)as.character(r[[j]]),character(1))
  }),stringsAsFactors=FALSE)
  names(frame)=gsub('"','',sub('.* AS ','',fields),fixed=TRUE)
  replacement=dbreg:::dbglm_values_sql(native,frame,format='values')
  statement=paste0(substr(statement,1,m_pos-1),replacement,substr(sub_stmt,end_pos+1,nchar(sub_stmt)))
 }
 statement=gsub('monotonically_increasing_id()','ROW_NUMBER() OVER ()',statement,fixed=TRUE)
 statement=gsub(' EXCEPT',' EXCLUDE',statement,fixed=TRUE)
 sub('^DROP TEMP TABLE','DROP TABLE',statement)
}
methods::setMethod('dbGetQuery',c('dbglm_sql_fixture','character'),function(conn,statement,...)dbGetQuery(conn@native,fixture_sql(statement),...))
methods::setMethod('dbExecute',c('dbglm_sql_fixture','character'),function(conn,statement,...)dbExecute(conn@native,fixture_sql(statement),...))
native=dbConnect(duckdb::duckdb(shared_home=FALSE),config=list(temp_directory='',threads='2'))
remote=methods::new('dbglm_sql_fixture',native=native)
r=dbreg:::dbglm_get_query(remote,"SELECT isfinite(1.0) AS f,isfinite('Infinity'::DOUBLE) AS i,isfinite('NaN'::DOUBLE) AS n,FSUM(x) AS s FROM (VALUES(1.0),(2.0)) t(x)")
expect_identical(as.logical(r[1,1:3]),c(TRUE,FALSE,FALSE))
expect_equal(r$s,3)
# Quoted identifiers and SQL string literals survive dialect translation.
r=dbreg:::dbglm_get_query(remote,"SELECT \"FSUM\",'FSUM isfinite EXCLUDE VARCHAR' AS text FROM (SELECT 7 AS \"FSUM\") t")
expect_equal(r$FSUM,7)
expect_identical(r$text,'FSUM isfinite EXCLUDE VARCHAR')
v=data.frame(k=c("a'b",'FSUM'),x=c(1e-25,-1e25))
s=dbreg:::dbglm_values_sql(remote,v)
expect_equal(dbGetQuery(remote,s),v,check.attributes=FALSE)
# Extreme doubles and escaped keys survive JSON and the SQL parser exactly.
edge=data.frame(k=c("a'b",'x\\y','\n',intToUtf8(7),'é','東京'),
 x=c(.65313212622670436,.Machine$double.xmax,.Machine$double.xmin,1e-300,-1e300,.1))
expect_identical(dbGetQuery(remote,dbreg:::dbglm_values_sql(remote,edge))$k,edge$k)
expect_identical(dbGetQuery(remote,dbreg:::dbglm_values_sql(remote,edge))$x,edge$x)
expect_equal(dbGetQuery(remote,dbreg:::dbglm_values_sql(remote,v,format='values')),v,check.attributes=FALSE)
expect_null(dbreg:::dbglm_values_sql(remote,v,max_bytes=1))
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
