library(dbreg)
library(DBI)
library(fixest)
source(system.file("tinytest/helper_glm_reference.R",package="dbreg"),local=TRUE)

(function() {
  con = dbConnect(duckdb::duckdb(),
    config=list(temp_directory="",memory_limit="2GiB",threads="2"))
  on.exit(dbDisconnect(con,shutdown=TRUE))
  eta = sort(unique(c(seq(-12,12,length.out=2401),-1000,1000,
    -sqrt(50),sqrt(50),-sqrt(50)+1e-9,sqrt(50)-1e-9)))
  dbWriteTable(con,"kernel_grid",data.frame(eta=eta))
  k = dbreg:::dbglm_family_sql("probit")
  got = dbGetQuery(con,paste("SELECT",k$mu,"AS mu,",k$derivative,"AS d FROM kernel_grid"))
  expect_true(max(abs(got$mu-binomial("probit")$linkinv(eta)))<8e-16)
  expect_true(max(abs(got$d-binomial("probit")$mu.eta(eta)))<8e-16)

  set.seed(20261023)
  n=1600L
  dat=data.frame(x=rnorm(n),z=rnorm(n),g=sample(1:17,n,TRUE),h=sample(1:13,n,TRUE),
    cluster=sample(1:31,n,TRUE),freq=sample(1:3,n,TRUE),off=runif(n,-.15,.15))
  families=list(logit=binomial("logit"),probit=binomial("probit"),
    cloglog=binomial("cloglog"),poisson=poisson())
  for(family in families) {
    predictor=-.25+.3*dat$x-.2*dat$z+dat$off+.12*sin(dat$g)+.15*cos(dat$h)
    dat$y=if(family$family=="poisson")rpois(n,family$linkinv(predictor)) else
      rbinom(n,1,family$linkinv(predictor))
    dbWriteTable(con,"family_fixture",dat,overwrite=TRUE)
    plain=dbglm(y~x+z,conn=con,table="family_fixture",family=family,offset="off",
      weights="freq",ssc=list(K.adj=FALSE))
    ref=glm(y~x+z,data=dat,family=family,offset=off,weights=freq,
      control=glm.control(epsilon=1e-13,maxit=100))
    expect_true(max(abs(coef(plain)-coef(ref)))<1e-6)
    expect_true(max(abs(vcov(plain)-vcov(ref)))<1e-6)
    expect_true(abs(as.numeric(logLik(plain)-logLik(ref)))<1e-6)
    expect_true(abs(plain$deviance-deviance(ref))<1e-6)
    expect_true(max(abs(predict(plain,dat,type="response")-predict(ref,dat,type="response")))<1e-6)
    expect_identical(generics::glance(plain)$family,family$family)
    for(covariance in list("iid","hc1",~cluster)) {
      got=dbglm(y~x+z|g+h,conn=con,table="family_fixture",family=family,
        weights="freq",offset="off",vcov=covariance,fixef_rm="none")
      reference=dbglm_reference_feglm(y~x+z|g+h,data=dat,family=family,weights=~freq,
        offset=~off,vcov=if(identical(covariance,"hc1"))"hetero" else covariance,fixef.rm="none")
      order=match(names(coef(got)),names(coef(reference)))
      expect_true(!anyNA(order))
      expect_true(max(abs(coef(got)-coef(reference)[order]))<1e-6)
      expect_true(max(abs(vcov(got)-vcov(reference)[order,order]))<1e-6)
      expect_true(abs(as.numeric(logLik(got)-logLik(reference)))/n<1e-9)
      expect_true(got$converged && got$score_max<=got$control$score_tol)
      expect_identical(got$nobs_orig,as.numeric(n))
    }
  }

  # A positive constant group is a finite Poisson intercept. Zero groups have
  # infinite coefficients and are removed under the usual FE pruning policy.
  pd=data.frame(x=rnorm(90),g=rep(1:3,each=30))
  pd$y=c(rep(c(0,1,2),10),rep(2,30),rep(0,30))
  fit=dbglm(y~x|g,data=pd,family=poisson())
  ref=dbglm_reference_feglm(y~x|g,data=pd,family=poisson())
  expect_equal(nobs(fit),60)
  expect_true(max(abs(coef(fit)-coef(ref)))<1e-6)
  expect_true(max(abs(vcov(fit)-vcov(ref)))<1e-6)
  expect_error(dbglm(y~x|g,data=pd,family=poisson(),fixef_rm="none"),"separation")
  constant=dbglm(y~1,data=data.frame(y=rep(3,30)),family=poisson(),ssc=list(K.adj=FALSE))
  expect_equal(unname(coef(constant)),log(3),tolerance=1e-10)
  expect_error(dbglm(y~1,data=data.frame(y=rep(0,10)),family=poisson()),"finite")
  expect_error(dbglm(y~x,data=transform(pd,y=-1),family=poisson()),"nonnegative integers")
  expect_error(dbglm(y~x,data=transform(pd,y=.5),family=poisson()),"nonnegative integers")
  expect_error(dbglm(y~x,data=pd,family=poisson(),separation="prediction"),"binomial")
  expect_error(dbglm(y~x,data=pd,family=poisson("identity")),"family")
  expect_error(dbglm(y~x,data=pd,family=binomial("cauchit")),"family")
  expect_error(dbglm_control(projection_method="unsupported"),"arg")
  expect_false(any(startsWith(dbListTables(con),"dbglm_")))
})()
