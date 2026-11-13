local({
  set.seed(11071)
  n=144L
  d=data.frame(g=rep(1:3,48),f=rep(c("a","b","c"),48),x=rnorm(n),w=runif(n,.4,2),
    cl=rep(1:12,12),fe=rep(1:8,18))
  d[["y"]]=.7+.5*d[["x"]]+.6*(d[["g"]]==3)+sin(d[["fe"]])+rnorm(n)
  ref=d;ref[["g"]]=relevel(factor(ref[["g"]]),ref="2");ref[["f"]]=relevel(factor(ref[["f"]]),ref="b")
  for(formula in list(y~g+x,y~g:x+x,y~0+g+x,y~f+x)) {
    source=if("g"%in%all.vars(formula))"g"else"f"
    refs=setNames(list(if(source=="g")"2"else"b"),source)
    for(kind in list("iid","hc1",~cl)) for(strategy in c("moments","compress")) {
      oracle=fixest::feols(formula,ref,weights=~w,vcov=kind,notes=FALSE)
      fit=dbreg(formula,data=d,weights="w",factor_refs=refs,vcov=kind,strategy=strategy)
      expect_equal(coef(fit),coef(oracle)[names(coef(fit))],tolerance=1e-8)
      expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(oracle)[names(coef(fit)),names(coef(fit)),drop=FALSE]),tolerance=1e-8)
      expect_equal(predict(fit,d),as.numeric(predict(oracle,ref)),tolerance=1e-8)
      if(strategy=="compress")expect_equal(sort(predict(fit)),sort(as.numeric(predict(oracle,ref))[!duplicated(model.matrix(oracle))]),tolerance=1e-8)
    }
  }
  for(strategy in c("demean","compress")) {
    oracle=fixest::feols(y~g+x|fe,ref,weights=~w,vcov=~cl,fixef.rm="none",ssc=fixest::ssc(K.fixef="full",K.exact=TRUE))
    fit=dbreg(y~g+x|fe,data=d,weights="w",factor_refs=c(g="2"),vcov=~cl,strategy=strategy)
    expect_equal(coef(fit),coef(oracle)[names(coef(fit))],tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(oracle)[names(coef(fit)),names(coef(fit)),drop=FALSE]),tolerance=1e-8)
  }
  # References persist independently of newdata factor order or contrast options.
  fit=dbreg(y~f+x,data=d,factor_refs=c(f="b"),strategy="compress")
  changed=d;changed[["f"]]=factor(changed[["f"]],levels=c("c","a","b"))
  expect_equal(predict(fit,changed),predict(fit,d),tolerance=1e-12)
  saved=unserialize(serialize(fit,NULL));expect_equal(predict(saved,d),predict(fit,d),tolerance=1e-12)
  bad=d[1:3,];bad[["f"]]=c("unknown",NA,"c")
  expect_equal(predict(saved,bad),c(NA_real_,NA_real_,predict(saved,d[3,,drop=FALSE])),tolerance=1e-12)
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  for(refs in list(list(g=NA),list(g=c(1,2)),list(unknown=1),list(2),setNames(list(1,2),c("g","g")),list(g=99))) {
    expect_error(dbreg(y~g+x,conn=conn,data=d,factor_refs=refs),"factor_refs|Reference level")
    expect_identical(DBI::dbListTables(conn),character())
  }
  # A reference on the response is invalid in the standalone compiler too.
  DBI::dbWriteTable(conn,"factor_fixture",d)
  expect_error(sql_model_matrix(y~g+x,conn,"factor_fixture",factor_refs=c(y=1)),"factor_refs")
  design=sql_model_matrix(~g+x,conn,"factor_fixture",factor_refs=c(g="2"))
  actual=DBI::dbGetQuery(conn,paste("SELECT",paste(design[["select_exprs"]],collapse=","),"FROM factor_fixture"))
  expect_equal(unname(as.matrix(actual)),unname(model.matrix(~g+x,ref)[,-1,drop=FALSE]),tolerance=1e-12)
  # Unsafe source names are literals in both factor references and i() keys.
  unsafe=d;names(unsafe)[names(unsafe)=="g"]="group name"
  DBI::dbWriteTable(conn,"unsafe_fixture",unsafe)
  for(formula in list(~`group name`,~i(`group name`,ref=2))) {
    design=sql_model_matrix(formula,conn,"unsafe_fixture",factor_refs=if(length(formula[[2]])==1L)setNames(list(2),"group name")else NULL)
    expect_equal(length(design[["aliases"]]),2L)
    values=DBI::dbGetQuery(conn,paste("SELECT",paste(design[["select_exprs"]],collapse=","),"FROM unsafe_fixture"))
    expect_equal(unname(as.matrix(values)),cbind(as.numeric(d[["g"]]==1),as.numeric(d[["g"]]==3)),tolerance=1e-12)
  }
})

local({
  # Numeric keys are matched before printing, preserving adjacent doubles.
  set.seed(11072)
  keys=1+c(0,1,2)*.Machine[["double.eps"]];index=rep(1:3,32)
  d=data.frame(g=keys[index],x=rnorm(96),y=rnorm(96))
  X=cbind(1,as.numeric(index==1),as.numeric(index==3),d[["x"]]);beta=as.numeric(qr.solve(X,d[["y"]]))
  for(strategy in c("moments","compress")) {
    fit=dbreg(y~g+x,data=d,factor_refs=list(g=keys[2]),strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    expect_false(anyDuplicated(names(coef(fit)))>0L)
    newer=d[1:3,];newer[["g"]][1]=1+3*.Machine[["double.eps"]]
    expect_true(is.na(predict(fit,newer)[1]))
  }
})

local({
  set.seed(11073)
  d = data.frame(g = rep(1 + .Machine[["double.eps"]], 40), x = rnorm(40), y = rnorm(40))
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y~g+x, data=d, factor_refs=list(g=d[["g"]][1]), strategy=strategy)
    expect_equal(predict(fit,d), as.numeric(predict(lm(y~x,d),d)), tolerance=1e-8)
    unknown = d[1:3,]; unknown[["g"]][1]=1; unknown[["g"]][2]=NA_real_
    expect_equal(predict(fit,unknown), c(NA_real_,NA_real_,predict(fit,d[3,,drop=FALSE])), tolerance=1e-8)
  }
})

local({
  set.seed(11074)
  n=120L; index=rep(1:3,40)
  d=data.frame(g=index,x=rnorm(n),y=rnorm(n),fe=rep(1:8,15))
  ref=d;ref[["g"]]=relevel(factor(ref[["g"]]),ref="2")
  fit=dbreg(y~g+x|fe,data=d,factor_refs=c(g=2),strategy="mundlak")
  oracle=fixest::feols(y~g+x|fe,ref,fixef.rm="none")
  expect_equal(coef(fit),coef(oracle)[names(coef(fit))],tolerance=1e-8)
  for(kind in c("logical","date")) {
    values=if(kind=="logical")rep(c(FALSE,TRUE),60)else as.Date("2020-01-01")+index-1
    d[["g"]]=values
    reference=if(kind=="logical")TRUE else as.Date("2020-01-02")
    selected=if(kind=="logical")FALSE else as.Date(c("2020-01-01","2020-01-03"))
    X=do.call(cbind,c(list(rep(1,n)),lapply(selected,function(key)as.numeric(values==key)),list(d[["x"]])))
    beta=as.numeric(qr.solve(X,d[["y"]]))
    for(strategy in c("moments","compress")) {
      fit=dbreg(y~g+x,data=d,factor_refs=list(g=reference),strategy=strategy)
      expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
      expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    }
  }
  conn=DBI::dbConnect(duckdb::duckdb(bigint="numeric",config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  keys=c("-9223372036854775808","9007199254740993","9223372036854775807")
  d[["g"]]=keys[index]
  DBI::dbWriteTable(conn,"factor_text",d)
  DBI::dbExecute(conn,"CREATE TABLE signed_factor AS SELECT CAST(g AS BIGINT) AS g,x,y FROM factor_text")
  X=cbind(1,as.numeric(index==1),as.numeric(index==3),d[["x"]]);beta=as.numeric(qr.solve(X,d[["y"]]))
  for(strategy in c("moments","compress")) {
    fit=dbreg(y~g+x,conn=conn,table="signed_factor",factor_refs=list(g=keys[2]),strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
  }
})
