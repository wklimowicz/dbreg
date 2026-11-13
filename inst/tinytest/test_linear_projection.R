local({
  set.seed(11091); n=180L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=sample(1:12,n,TRUE),h=sample(1:9,n,TRUE),k=sample(1:5,n,TRUE),w=runif(n,.2,2),cl=rep(1:15,12))
  d[["y"]]=.7+.5*d[["x"]]-.3*d[["z"]]+sin(d[["g"]])+rnorm(n)
  meanw=function(x,g,w)ave(seq_along(x),g,FUN=function(i)weighted.mean(x[i],w[i]))
  for(scale in c(1,1e-200,1e200)) {
    data=d;data[["w"]]=data[["w"]]*scale
    X=cbind(1,d[["x"]],d[["z"]],meanw(d[["x"]],d[["g"]],d[["w"]]),meanw(d[["z"]],d[["g"]],d[["w"]]))
    fit=dbreg(y~x+z|g,data=data,weights="w",strategy="mundlak",vcov=~cl)
    beta=as.numeric(qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]])))
    expect_equal(as.numeric(coef(fit,fe=TRUE)),beta,tolerance=1e-8)
    expect_equal(predict(fit,data),as.numeric(X%*%beta),tolerance=1e-8)
  }
  for(formula in list(y~x+z|g,y~x+z|g+h,y~x+z|g+h+k)) {
    fit=dbreg(formula,data=d,weights="w",strategy="demean")
    rhs=if(length(all.vars(formula))==4L)~x+z+factor(g)else if(length(all.vars(formula))==5L)~x+z+factor(g)+factor(h)else ~x+z+factor(g)+factor(h)+factor(k)
    X=model.matrix(rhs,d)
    beta=qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]]))
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-7)
    no_y=d;no_y[["y"]]=NULL
    Z=X[,setdiff(seq_len(ncol(X)),match(c("x","z"),colnames(X))),drop=FALSE]
    U=as.matrix(d[c("x","z")]);within=U-Z%*%qr.solve(Z*sqrt(d[["w"]]),U*sqrt(d[["w"]]))
    expect_equal(suppressMessages(predict(fit,no_y)),as.numeric(within%*%coef(fit)),tolerance=1e-7)
    changed=d[rev(1:30),];changed[["x"]][2]=NA;changed[["w"]][3]=NA
    found=predict(fit,changed)
    expect_equal(length(found),30L);expect_true(all(is.na(found[2:3])))
    expect_equal(found[-c(2,3)],predict(fit,changed[-c(2,3),]),tolerance=1e-9)
  }
  expect_error(predict(fit,transform(d,w=-1)),"non-negative")
  empty=d[FALSE,];expect_identical(predict(fit,empty),numeric())
  expect_error(dbreg:::demean_prediction_data(as.matrix(d[c("x","z")]),list(d[["g"]],d[["h"]]),d[["w"]],tol=1e-14,max_iter=1L),"did not converge")
})

local({
  set.seed(11092);n=160L
  d=data.frame(x=rnorm(n),z=rnorm(n),w=runif(n,.2,2),g=rep(1:8,20),cl=rep(1:10,16))
  d[["y"]]=.5*d[["x"]]+sin(d[["g"]])+rnorm(n)
  means=function(x)ave(seq_len(n),d[["g"]],FUN=function(i)weighted.mean(x[i],d[["w"]][i]))
  X=cbind(1,d[["x"]],d[["z"]],means(d[["x"]]),means(d[["z"]]))
  beta=as.numeric(qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]])))
  for(name in c("weights","dbreg_bar1_1")) {
    data=d;names(data)[names(data)=="x"]=name
    formula=as.formula(paste("y ~",paste(c(name,"z"),collapse=" + "),"| g"))
    fit=dbreg(formula,data=data,weights="w",strategy="mundlak")
    expect_equal(as.numeric(coef(fit,fe=TRUE)),beta,tolerance=1e-8)
    expect_equal(predict(fit,data),as.numeric(X%*%beta),tolerance=1e-8)
  }
})

local({
  set.seed(11093);n=150L
  d=data.frame(x=rnorm(n),g=sample(1:10,n,TRUE),f=rep(1:3,50),w=runif(n,.2,2),cl=rep(1:15,10))
  d[["y"]]=.7*d[["x"]]+.2*d[["x"]]^2+sin(d[["g"]])+rnorm(n)
  raw=model.matrix(~x+I(x^2)+factor(f),d)[,-1,drop=FALSE]
  bars=apply(raw,2,function(x)ave(seq_len(n),d[["g"]],FUN=function(i)weighted.mean(x[i],d[["w"]][i])))
  X=cbind(1,raw,bars);w=d[["w"]];beta=as.numeric(qr.solve(X*sqrt(w),d[["y"]]*sqrt(w)))
  e=d[["y"]]-as.numeric(X%*%beta);bread=solve(crossprod(X,w*X))
  for(kind in list("iid","hc1",~cl)) {
    meat=if(identical(kind,"hc1"))crossprod(X*(w*e))else if(inherits(kind,"formula"))crossprod(rowsum(X*(w*e),d[["cl"]]))else NULL
    V=if(is.null(meat))bread*sum(w*e^2)/(n-ncol(X))else bread%*%meat%*%bread*if(inherits(kind,"formula"))15/14*149/(n-ncol(X))else n/(n-ncol(X))
    fit=dbreg(y~x+I(x^2)+factor(f)|g,data=d,weights="w",vcov=kind,strategy="mundlak")
    expect_equal(as.numeric(coef(fit,fe=TRUE)),beta,tolerance=1e-8)
    expect_equal(as.numeric(vcov(fit,fe=TRUE)),as.numeric(V),tolerance=1e-8)
    expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-8)
    bad=d[1:10,];bad[["x"]][2]=NA;bad[["w"]][3]=NA
    found=predict(fit,bad);expect_true(all(is.na(found[2:3])))
    expect_equal(found[-c(2,3)],predict(fit,bad[-c(2,3),]),tolerance=1e-10)
    missing=bad;missing[["w"]]=NULL;expect_error(predict(fit,missing),"missing required")
  }
})

local({
  set.seed(11094);n=144L;group=rep(1:12,12)
  d=data.frame(x=rnorm(n),z=rnorm(n),w=runif(n,.2,2),g=1+group*.Machine[["double.eps"]])
  d[["y"]]=.5*d[["x"]]-.3*d[["z"]]+sin(group)+rnorm(n)
  Z=model.matrix(~factor(group));X=cbind(Z,as.matrix(d[c("x","z")]))
  beta=qr.solve(X*sqrt(d[["w"]]),d[["y"]]*sqrt(d[["w"]]))
  fit=dbreg(y~x+z|g,data=d,weights="w",strategy="demean")
  expect_equal(predict(fit,d),as.numeric(X%*%beta),tolerance=1e-7)
  dropped=d;dropped[["w"]][group==1]=0
  found=predict(fit,dropped)
  expect_true(all(is.na(found[group==1])))
  expect_equal(found[group!=1],predict(fit,dropped[group!=1,]),tolerance=1e-9)
  all_zero=d;all_zero[["w"]]=0;expect_true(all(is.na(predict(fit,all_zero))))
  missing=d;missing[["g"]][1]=NA;expect_true(is.na(predict(fit,missing)[1]))
})

local({
  set.seed(11096);n=180L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=sample(1:10,n,TRUE),h=sample(1:8,n,TRUE),w=runif(n,.2,2))
  d[["endog"]]=.7*d[["z"]]+rnorm(n);d[["y"]]=.5*d[["x"]]+.8*d[["endog"]]+sin(d[["g"]])+rnorm(n)
  F=model.matrix(~factor(g)+factor(h),d)
  U=as.matrix(d[c("x","endog","z","y")]);w=d[["w"]]
  projected=U-F%*%qr.solve(F*sqrt(w),U*sqrt(w))
  X=projected[,1:2];Z=projected[,c(1,3)];y=projected[,4]
  P=crossprod(X,w*Z)%*%solve(crossprod(Z,w*Z));beta=as.numeric(solve(P%*%crossprod(Z,w*X),P%*%crossprod(Z,w*y)))
  fit=dbivreg(y~x|g+h|endog~z,data=d,weights="w",strategy="demean")
  expect_equal(as.numeric(coef(fit)),beta,tolerance=1e-8)
  expect_equal(predict(fit,d),as.numeric(X%*%beta+d[["y"]]-y),tolerance=1e-7)
  no_y=d;no_y[["y"]]=NULL
  expect_equal(suppressMessages(predict(fit,no_y)),as.numeric(X%*%beta),tolerance=1e-7)
})
