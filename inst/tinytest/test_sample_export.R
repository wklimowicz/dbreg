library(DBI)
local({
  conn=DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(11110);n=160L
  d=data.frame(id=seq_len(n),x=as.numeric(rep(c(-2,-1,1,2),40)),z=as.numeric(rep(0:3,40)),g=as.numeric(rep(1:8,each=20)),w=runif(n,.2,2),cl=as.numeric(rep(1:10,16)),o=rnorm(n),pair=rep(c("a|b",""),80),day=as.Date('2026-01-01')+seq_len(n),flag=rep(c(TRUE,FALSE),80))
  d[["endog"]]=.7*d[["z"]]+rnorm(n);d[["y"]]=d[["o"]]+.5*d[["x"]]+.7*d[["endog"]]+rnorm(n)
  DBI::dbWriteTable(conn,"sample_raw",d)
  DBI::dbExecute(conn,'CREATE TEMP TABLE sample_source AS SELECT *, CAST(\'9007199254740000\' AS BIGINT)+CAST(1100+CEIL(id/2.0) AS BIGINT) AS "row key", pair AS "dbreg_key_count__" FROM sample_raw')
  DBI::dbExecute(conn,"UPDATE sample_source SET x=CAST('NaN' AS DOUBLE) WHERE id=2")
  DBI::dbExecute(conn,"UPDATE sample_source SET y=NULL WHERE id=3")
  DBI::dbExecute(conn,"UPDATE sample_source SET w=0 WHERE id=4")
  DBI::dbExecute(conn,"UPDATE sample_source SET cl=NULL WHERE id=5")
  DBI::dbExecute(conn,"UPDATE sample_source SET o=CAST('Infinity' AS DOUBLE) WHERE id=6")
  keys=c("row key","dbreg_key_count__","day","flag")
  before=sort(DBI::dbListTables(conn))
  keep=d[["x"]]>-2 & !d[["id"]] %in% 2:6
  expected=data.frame(key=paste0("900719925474",1100+ceiling(d[["id"]][keep]/2)),pair=d[["pair"]][keep],day=d[["day"]][keep],flag=d[["flag"]][keep])
  expected=expected[order(expected[["key"]],expected[["pair"]]),];rownames(expected)=NULL
  read_export=function(){
    values=DBI::dbGetQuery(conn,'SELECT CAST("row key" AS VARCHAR) AS key,"dbreg_key_count__" AS pair,day,flag FROM "sample keys" ORDER BY "row key","dbreg_key_count__"')
    rownames(values)=NULL;values
  }
  for(estimator in c("ols","iv"))for(strategy in if(estimator=="ols")c("moments","compress","demean","mundlak")else c("moments","compress","demean")){
    has_fe=strategy %in% c("compress","demean","mundlak") && estimator=="ols" || strategy=="demean"
    formula=if(estimator=="ols"){if(has_fe)y~x+z|g else y~x+z}else{if(has_fe)y~x|g|endog~z else y~x|endog~z}
    fun=if(estimator=="ols")dbreg else dbivreg
    args=list(fml=formula,conn=conn,table="sample_source",weights=~w,offset=~o,subset=~x > -2,vcov=~cl,strategy=strategy)
    fit=do.call(fun,c(args,list(row_id=keys,sample_table="sample keys",diagnostics=TRUE)))
    plain=do.call(fun,args)
    expect_equal(coef(fit,fe=TRUE),coef(plain,fe=TRUE),tolerance=1e-9)
    expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(vcov(plain,fe=TRUE)),tolerance=1e-9)
    expect_identical(fit[["row_id"]],keys);expect_identical(fit[["sample_table"]],"sample keys")
    expect_equal(nobs(fit),sum(keep));expect_equal(fit[["nobs_input"]],n)
    expect_identical(read_export(),expected)
    schema=DBI::dbGetQuery(conn,'DESCRIBE "sample keys"')
    expect_identical(schema[["column_type"]],c("BIGINT","VARCHAR","DATE","BOOLEAN"))
    expect_equal(sort(DBI::dbListTables(conn)),sort(c(before,"sample keys")))
    expect_error(do.call(fun,c(args,list(row_id=keys,sample_table="sample keys"))),"already exists")
    expect_identical(read_export(),expected)
    DBI::dbExecute(conn,'DROP TABLE "sample keys"')
  }
  # Invalid keys are rejected on the full source, including excluded rows.
  expect_error(dbreg(y~x,conn=conn,table="sample_source",row_id="row key",sample_table="sample keys",subset=~x > -2),"unique and nonmissing")
  expect_equal(sort(DBI::dbListTables(conn)),before)
  DBI::dbExecute(conn,'CREATE TEMP TABLE missing_key AS SELECT * FROM sample_source')
  DBI::dbExecute(conn,'UPDATE missing_key SET "row key"=NULL WHERE id=2')
  expect_error(dbreg(y~x,conn=conn,table="missing_key",row_id=keys,sample_table="sample keys",subset=~id > 10),"unique and nonmissing")
  DBI::dbExecute(conn,'DROP TABLE missing_key')
  # A DOUBLE key exercises true SQL numeric missingness.
  DBI::dbExecute(conn,"CREATE TEMP TABLE missing_key AS SELECT CASE WHEN id=2 THEN CAST('NaN' AS DOUBLE) ELSE CAST(id AS DOUBLE) END AS bad_id,* FROM sample_source")
  expect_error(dbreg(y~x,conn=conn,table="missing_key",row_id="bad_id",sample_table="sample keys",subset=~id > 10),"unique and nonmissing")
  DBI::dbExecute(conn,'DROP TABLE missing_key')
  expect_error(dbreg(y~x,conn=conn,data=d,row_id=c("id","id")),"unique source key columns")
  expect_error(dbreg(y~x,conn=conn,data=d,sample_table="sample keys"),"row_id")
  expect_error(dbreg(y~x,data=d,row_id="id",sample_table="sample keys"),"caller-owned connection")
  expect_error(dbivreg(y~x|endog~z,conn=conn,data=d,row_id="id",sample_table="sample keys",sql_only=TRUE),"sql_only or data_only")
  expect_error(dbreg(y~x,conn=conn,data=d,row_id="id",sample_table="sample keys",data_only=TRUE),"sql_only or data_only")
  expect_equal(sort(DBI::dbListTables(conn)),before)
  # A caller-provided data.frame is registered only for the duration of the fit.
  fit=dbreg(y~x,conn=conn,data=d,row_id="id",sample_table="sample keys")
  expect_equal(DBI::dbGetQuery(conn,'SELECT COUNT(*) AS n FROM "sample keys"')[["n"]],n)
  expect_equal(sort(DBI::dbListTables(conn)),sort(c(before,"sample keys")))
  DBI::dbExecute(conn,'DROP TABLE "sample keys"')
  expect_true(DBI::dbIsValid(conn))
  # Signed BIGINT endpoints and adjacent IEEE doubles must not be rounded or
  # concatenated when deciding whether a composite key is unique.
  DBI::dbExecute(conn, paste0("CREATE TEMP TABLE exact_source AS SELECT *, ",
    "CASE WHEN id%2=0 THEN CAST('-9223372036854775808' AS BIGINT) ELSE CAST('9223372036854775807' AS BIGINT) END AS big, ",
    "CASE WHEN id%2=0 THEN 1.0::DOUBLE ELSE 1.0000000000000002::DOUBLE END AS adjacent FROM sample_raw"))
  fit=dbreg(y~x,conn=conn,table="exact_source",row_id=c("big","adjacent","day"),sample_table="sample keys")
  expect_equal(nobs(fit),n)
  expect_equal(DBI::dbGetQuery(conn,'SELECT COUNT(DISTINCT adjacent) AS n FROM "sample keys"')[["n"]],2)
  expect_identical(sort(DBI::dbGetQuery(conn,'SELECT DISTINCT CAST(big AS VARCHAR) AS key FROM "sample keys"')[["key"]]),c("-9223372036854775808","9223372036854775807"))
  DBI::dbExecute(conn,'DROP TABLE "sample keys"')
  DBI::dbExecute(conn,"CREATE TEMP TABLE unique_big AS SELECT *, CAST('9007199254740992' AS BIGINT)+CAST(id AS BIGINT) AS big FROM sample_raw")
  fit=dbreg(y~x,conn=conn,table="unique_big",row_id="big",sample_table="sample keys")
  expect_equal(nobs(fit),n)
  expect_identical(DBI::dbGetQuery(conn,'SELECT CAST(big AS VARCHAR) AS key FROM "sample keys" ORDER BY big')[["key"]],paste0("900719925474",sprintf("%04d",992+seq_len(n))))
  DBI::dbExecute(conn,'DROP TABLE "sample keys"')
  # Each of these two-row sources has unique keys that lossy conversion can
  # collapse. Use an intercept-only fit so there is no artificial extra key.
  for(source in c(
    "(SELECT * FROM (VALUES (1.0::DOUBLE,1.0),(1.0000000000000002::DOUBLE,2.0)) AS v(k,y)) exact",
    "(SELECT * FROM (VALUES ('a|b','c',1.0),('a','b|c',2.0)) AS v(a,b,y)) exact")){
    keys=if(grepl("v(k,y)",source,fixed=TRUE))"k" else c("a","b")
    fit=dbreg(y~1,conn=conn,table=source,row_id=keys,sample_table="sample keys")
    expect_equal(as.numeric(coef(fit)),1.5)
    expect_equal(DBI::dbGetQuery(conn,'SELECT COUNT(*) AS n FROM "sample keys"')[["n"]],2)
    DBI::dbExecute(conn,'DROP TABLE "sample keys"')
  }
  # Lazy sources infer the caller session; SQL pattern characters in names
  # must be literal when checking for an existing export.
  if(requireNamespace("dplyr",quietly=TRUE)&&requireNamespace("dbplyr",quietly=TRUE)){
    lazy=dplyr::tbl(conn,"sample_raw")
    fit=dbreg(y~x,data=lazy,row_id="id",sample_table="export_%| name")
    expect_equal(DBI::dbGetQuery(conn,'SELECT COUNT(*) AS n FROM "export_%| name"')[["n"]],n)
    expect_error(dbreg(y~x,data=lazy,row_id="id",sample_table="EXPORT_%| NAME"),"already exists")
    expect_true(DBI::dbIsValid(conn))
    DBI::dbExecute(conn,'DROP TABLE "export_%| name"')
  }
})
