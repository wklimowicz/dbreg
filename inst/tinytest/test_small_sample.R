library(DBI)
local({
  set.seed(11114)
  n=180L
  d=data.frame(x=rnorm(n),z=rnorm(n),g=rep(1:9,each=20),cl=rep(1:6,30),w=runif(n,.2,2))
  d[["endog"]]=.8*d[["z"]]+rnorm(n)
  d[["y"]]=.3*d[["x"]]+.6*d[["endog"]]+sin(d[["g"]])+rnorm(n)
  expected_policy=list(K.adj=TRUE,K.fixef="full",K.exact=FALSE,G.adj=TRUE,G.df="min",t.df="min")
  expect_identical(dbreg:::normalize_ssc("full"),expected_policy)
  expected_policy[["K.fixef"]]="nonnested"
  expect_identical(dbreg:::normalize_ssc("nested"),expected_policy)
  expect_identical(dbreg:::normalize_ssc(NULL),expected_policy)
  expected_policy[["K.fixef"]]="none"
  expect_identical(dbreg:::normalize_ssc("none"),expected_policy)
  policies=list("full","nested","none",list(K.adj=FALSE,G.adj=FALSE),
    list(K.fixef="none",t.df=7),list(K.fixef="full",t.df="conventional"))
  for(estimator in c("ols","iv"))for(strategy in if(estimator=="ols")c("moments","compress","demean","mundlak")else c("moments","compress","demean")){
    has_fe=strategy %in% c("demean","mundlak") || estimator=="ols" && strategy=="compress"
    formula=if(estimator=="ols"){if(has_fe)y~x+z|g else y~x+z}else{if(has_fe)y~x|g|endog~z else y~x|endog~z}
    fun=if(estimator=="ols")dbreg else dbivreg
    for(kind in list("iid","hc1",~g,~cl))for(policy in policies){
      info=paste(estimator,strategy,paste(deparse(kind),collapse=""),paste(deparse(policy),collapse=""))
      normalized=dbreg:::normalize_ssc(policy)
      reference_policy=do.call(fixest::ssc,normalized)
      data=d;fml_ref=formula
      if(strategy=="mundlak"){
        meanw=function(x)ave(seq_len(nrow(d)),d[["g"]],FUN=function(i)weighted.mean(x[i],d[["w"]][i]))
        data[["mx"]]=meanw(d[["x"]]);data[["mz"]]=meanw(d[["z"]]);fml_ref=y~x+z+mx+mz
      }
      reference=fixest::feols(fml_ref,data,weights=~w,vcov=if(identical(kind,"hc1"))"hetero"else kind,
        ssc=reference_policy,fixef.rm="none",notes=FALSE)
      fit=fun(formula,data=d,weights=~w,strategy=strategy,vcov=kind,ssc=policy)
      full=strategy=="mundlak"
      ix=if(full)seq_along(coef(reference))else match(names(coef(fit)),sub("^fit_","",names(coef(reference))))
      expect_equal(as.numeric(vcov(fit,fe=full)),as.numeric(vcov(reference)[ix,ix,drop=FALSE]),tolerance=1e-7,info=info)
      ct=if(full)fit[["coeftable"]] else fit[["coeftable"]][names(coef(fit)),,drop=FALSE]
      expect_equal(as.numeric(ct[,"p.value"]),as.numeric(fixest::coeftable(reference)[ix,4]),tolerance=1e-7,info=info)
      expect_equal(attr(vcov(fit),"df.K"),fixest::degrees_freedom(reference,"k"),info=info)
      expect_equal(attr(vcov(fit),"df.t"),fixest::degrees_freedom(reference,"t"),info=info)
      expect_identical(attr(vcov(fit),"ssc"),normalized,info=info)
      if(!full)expect_equal(as.numeric(confint(fit)),as.numeric(as.matrix(confint(reference)[ix,,drop=FALSE])),tolerance=1e-7,info=info)
    }
  }
  # K.exact=FALSE uses a declared FE-count convention, including nested or
  # disconnected dummy systems; compression must not silently use another K.
  d[["h"]]=ceiling(d[["g"]]/3)
  for(strategy in c("compress","demean"))for(policy in policies)for(kind in list("iid",~h)){
    reference=fixest::feols(y~x|g+h,d,weights=~w,vcov=kind,
      ssc=do.call(fixest::ssc,dbreg:::normalize_ssc(policy)),fixef.rm="none",notes=FALSE)
    fit=dbreg(y~x|g+h,data=d,weights=~w,vcov=kind,strategy=strategy,ssc=policy)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(reference)),tolerance=1e-7)
    expect_equal(attr(vcov(fit),"df.K"),fixest::degrees_freedom(reference,"k"))
    expect_equal(as.numeric(fit[["coeftable"]][names(coef(fit)),"p.value"]),as.numeric(fixest::coeftable(reference)[,4]),tolerance=1e-7)
  }
  # Test df changes intervals, while the prediction noise still uses model df.
  fit=dbreg(y~x,data=d,vcov=~g,ssc=list(t.df=7),strategy="moments")
  X=model.matrix(~x,d);V=as.matrix(vcov(fit));point=as.numeric(X%*%coef(fit))
  se=sqrt(rowSums((X%*%V)*X));width=qt(.975,7)*se
  ci=predict(fit,d,interval="confidence")
  expect_equal(ci[["lwr"]],unname(point-width),tolerance=1e-9)
  expect_equal(ci[["upr"]],unname(point+width),tolerance=1e-9)
  pi=predict(fit,d,interval="prediction")
  noise=attr(vcov(fit),"rss")/fit[["df_residual"]]
  expect_equal(pi[["upr"]],unname(point+qt(.975,7)*sqrt(se^2+noise)),tolerance=1e-9)
  # A first-stage Wald test uses its own instrument-space correction and test df.
  iv=dbivreg(y~x|endog~z,data=d,vcov=~g,ssc=list(t.df=7),strategy="moments")
  first=dbreg(endog~x+z,data=d,vcov=~g,ssc=list(t.df=7),strategy="moments")
  expected=coef(first)[["z"]]^2/vcov(first)["z","z"]
  actual=iv[["diagnostics"]][["first_stage_wald"]][["endog"]]
  expect_equal(actual[["stat"]],expected,tolerance=1e-8)
  expect_equal(actual[["df2"]],7)
  expect_equal(actual[["p"]],pf(expected,1,7,lower.tail=FALSE),tolerance=1e-8)
  d[["z2"]]=rnorm(n)
  for(policy in policies)for(has_fe in c(FALSE,TRUE)){
    iv=dbivreg(if(has_fe)y~x|g|endog~z+z2 else y~x|endog~z+z2,data=d,weights=~w,vcov=~g,ssc=policy,strategy=if(has_fe)"demean"else"moments")
    first=fixest::feols(if(has_fe)endog~x+z+z2|g else endog~x+z+z2,d,weights=~w,vcov=~g,
      ssc=do.call(fixest::ssc,dbreg:::normalize_ssc(policy)),fixef.rm="none",notes=FALSE)
    b=coef(first)[c("z","z2")];V=vcov(first)[c("z","z2"),c("z","z2")]
    expected=as.numeric(crossprod(b,solve(V,b)))/2;df=fixest::degrees_freedom(first,"t")
    actual=iv[["diagnostics"]][["first_stage_wald"]][["endog"]]
    expect_equal(actual[["stat"]],expected,tolerance=1e-8)
    expect_equal(actual[["df2"]],df)
    expect_equal(actual[["p"]],pf(expected,2,df,lower.tail=FALSE),tolerance=1e-8)
  }
  for(policy in list("unknown",c("full","nested"),list(),list(foo=TRUE),list(K.exact=TRUE),list(K.adj=NA),list(t.df=Inf),list(t.df=0),list(G.df="bad"))){
    expect_error(dbreg(y~x,data=d,ssc=policy))
    expect_error(dbivreg(y~x|endog~z,data=d,ssc=policy))
  }
  expect_error(dbreg(y~x,data=transform(d,g=1),vcov=~g),"at least two clusters")
  expect_error(dbreg(y~x,data=d[1:2,],strategy="moments"),"residual degrees of freedom")
})
