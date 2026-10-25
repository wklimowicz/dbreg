library(dbreg)
library(fixest)
library(DBI)
ct <- dbConnect(duckdb::duckdb(config=list(temp_directory="")))
source(system.file('tinytest/helper_glm_reference.R',package='dbreg'),local=TRUE)
set.seed(271)
dt <- data.frame(x=rnorm(640),district=sample(1:8,640,TRUE),age=sample(1:5,640,TRUE),yr=sample(0:5,640,TRUE),a=sample(1:3,640,TRUE),o=runif(640,-.2,.2))
dt$yr2 <- dt$yr^2
dt$y <- rbinom(640,1,plogis(-.4+.25*dt$x-.01*dt$district*dt$yr))
for(link in c("cloglog","logit")) for(trend in c("yr","yr,yr2")) for(kind in c("iid","cluster")) {
 fm <- as.formula(paste0("y~x|district[",trend,"]+age"))
 vc <- if(kind=="cluster") ~district else "iid"
 m <- dbglm(fm,conn=ct,data=dt,family=binomial(link),weights="a",offset="o",vcov=vc)
 args <- list(fml=fm,data=dt,family=binomial(link),weights=~a,offset=~o,vcov=vc,glm.tol=3e-13,fixef.tol=1e-10,glm.iter=100)
 r <- do.call(dbglm_reference_feglm,args)
 expect_identical(names(coef(m)),names(coef(r)))
 expect_true(all(abs(coef(m)-coef(r))<=1e-7+1e-6*abs(coef(r))))
 expect_true(all(abs(sqrt(diag(vcov(m)))-sqrt(diag(vcov(r))))<=1e-8+1e-6*sqrt(diag(vcov(r)))))
 expect_true(all(abs(vcov(m)-vcov(r))<=1e-8+1e-5*abs(vcov(r))))
 expect_equal(as.numeric(logLik(m)),as.numeric(logLik(r)),tolerance=1e-8)
 expect_equal(m$n_params_model,r$nparams)
 expect_equal(m$n_params_vcov,attr(vcov(r,attr=TRUE),"df.K"))
 expect_true(m$score_max<=m$control$score_tol)
 expect_true(m$ap_residual<=m$control$fe_tol)
 expect_equal(length(dbListTables(ct)),0L)
}
# Constant trends in one district are harmless zero projection directions.
dt$yr[dt$district==1] <- 0
dt$yr2 <- dt$yr^2
fm <- y~x|district[yr,yr2]+age
m <- dbglm(family = stats::binomial("cloglog"), fm,conn=ct,data=dt)
r <- dbglm_reference_feglm(fm,dt,family=binomial('cloglog'))
expect_equal(coef(m),coef(r),tolerance=1e-7)
expect_equal(m$n_params_model,r$nparams)
expect_equal(length(dbListTables(ct)),0L)
# Source trend names must not collide with IRLS working columns.
dt$z <- dt$yr;dt$w <- dt$yr2
fm <- y~x|district[z,w]+age
m <- dbglm(family = stats::binomial("cloglog"), fm,conn=ct,data=dt)
r <- dbglm_reference_feglm(fm,dt,family=binomial('cloglog'))
expect_equal(coef(m),coef(r),tolerance=1e-7)
expect_equal(m$n_params_model,r$nparams)
expect_equal(length(dbListTables(ct)),0L)
dbDisconnect(ct,shutdown=TRUE)
