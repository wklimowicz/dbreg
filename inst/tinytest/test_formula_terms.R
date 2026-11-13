local({
  set.seed(11081);n=160L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=rep(1:4,40),f=rep(letters[1:4],each=40),w=runif(n,.4,2),cl=rep(1:16,10),fe=sample(rep(1:8,20)))
  d[["y"]]=.7+.5*d[["x"]]-.2*d[["x"]]^2+.1*d[["x"]]^3+sin(d[["g"]])+rnorm(n)
  cases=list(y~x+I(x^2)+I(x^3),y~factor(g)+x,y~I(x^2)*factor(g),y~I(x^2):f+f,
    y~0+factor(g)+I(x^2),y~x+factor(g):z,y~I(x^2)+I(z^3),y~factor(g)*f)
  for(formula in cases)for(kind in list("iid","hc1",~cl))for(strategy in c("moments","compress")) {
    oracle=lm(formula,d,weights=w)
    expected=coef(oracle)
    X=model.matrix(oracle);weights=d[["w"]];e=residuals(oracle);bread=solve(crossprod(X,weights*X));p=ncol(X)
    meat=if(identical(kind,"hc1"))crossprod(X*(weights*e))else if(inherits(kind,"formula"))crossprod(rowsum(X*(weights*e),d[["cl"]]))else NULL
    V=if(is.null(meat))bread*sum(weights*e^2)/(n-p)else bread%*%meat%*%bread*if(inherits(kind,"formula"))16/15*159/(n-p)else n/(n-p)
    dimnames(V)=rep(list(names(expected)),2L)
    fit=dbreg(formula,data=d,weights="w",vcov=kind,strategy=strategy)
    expect_equal(sort(names(coef(fit))),sort(names(expected)))
    expect_equal(coef(fit),expected[names(coef(fit))],tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(V[names(coef(fit)),names(coef(fit)),drop=FALSE]),tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(predict(oracle,d)),tolerance=1e-8)
    if(strategy=="compress")expect_equal(sort(predict(fit)),sort(as.numeric(predict(oracle,d))[!duplicated(model.matrix(oracle))]),tolerance=1e-8)
  }
  for(strategy in c("demean","compress")) {
    formula=y~x+I(x^2)+factor(g)|fe
    fit=dbreg(formula,data=d,weights="w",vcov=~cl,strategy=strategy)
    oracle=fixest::feols(formula,d,weights=~w,vcov=~cl,fixef.rm="none",ssc=fixest::ssc(K.fixef="full",K.exact=TRUE),notes=FALSE)
    expected=coef(oracle);names(expected)=sub("^I\\(I\\((.*)\\)\\)$","I(\\1)",names(expected))
    expect_equal(coef(fit),expected[names(coef(fit))],tolerance=1e-8)
  }
  # Literal source names resembling functions and generated aliases remain distinct.
  d[["I(x^2)"]]=rnorm(n);d[["factor(g)"]]=rnorm(n);d[["dbreg_power_1__"]]=rnorm(n)
  formula=y~I(x^2)+`I(x^2)`+factor(g)+`factor(g)`+dbreg_power_1__
  X=cbind(1,d[["x"]]^2,d[["I(x^2)"]],as.numeric(d[["g"]]==2),as.numeric(d[["g"]]==3),as.numeric(d[["g"]]==4),d[["factor(g)"]],d[["dbreg_power_1__"]])
  beta=as.numeric(qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]])))
  for(strategy in c("moments","compress")) {
    fit=dbreg(formula,data=d,weights="w",strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    saved=unserialize(serialize(fit,NULL));expect_equal(predict(saved,d),predict(fit,d),tolerance=1e-12)
  }
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="2")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE),add=TRUE)
  DBI::dbWriteTable(conn,"literal_sources",d)
  design=sql_model_matrix(formula,conn,"literal_sources")
  sql=paste("SELECT",paste(design[["select_exprs"]],collapse=","),"FROM literal_sources")
  expect_equal(unname(as.matrix(DBI::dbGetQuery(conn,sql))),unname(X[,-1,drop=FALSE]),tolerance=1e-12)
  for(formula in list(y~I(x^4),y~I(x+1),y~factor(g,levels=1:4),y~log(x))) {
    expect_error(dbreg(formula,data=d),"Only I|without additional|not supported")
  }
})

local({
  set.seed(11082);n=180L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=rep(1:3,60),w=runif(n,.4,2),cl=rep(1:15,12))
  d[["endog"]]=.7*d[["z"]]+.4*d[["z"]]^2+rnorm(n)
  d[["y"]]=.5*d[["x"]]-.3*d[["x"]]^2+.8*d[["endog"]]+sin(d[["g"]])+rnorm(n)
  X=model.matrix(~x+I(x^2)+factor(g)+endog,d);Z=model.matrix(~x+I(x^2)+factor(g)+z+I(z^2),d)
  w=d[["w"]];P=crossprod(X,w*Z)%*%solve(crossprod(Z,w*Z));bread=solve(P%*%crossprod(Z,w*X))
  beta=as.numeric(bread%*%P%*%crossprod(Z,w*d[["y"]]));e=d[["y"]]-as.numeric(X%*%beta)
  for(kind in list("iid","hc1",~cl))for(strategy in c("moments","compress")) {
    meat=if(identical(kind,"hc1"))crossprod(Z*(w*e))else if(inherits(kind,"formula"))crossprod(rowsum(Z*(w*e),d[["cl"]]))else NULL
    V=if(is.null(meat))bread*sum(w*e^2)/(n-ncol(X))else bread%*%P%*%meat%*%t(P)%*%bread*if(inherits(kind,"formula"))15/14*179/(n-ncol(X))else n/(n-ncol(X))
    fit=dbivreg(y~x+I(x^2)+factor(g)|endog~z+I(z^2),data=d,weights="w",vcov=kind,strategy=strategy)
    expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit)),as.numeric(V),tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    if(strategy=="compress")expect_equal(sort(predict(fit)),sort(as.numeric(X%*%beta)),tolerance=1e-8)
  }
})
