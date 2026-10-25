library(dbreg)
library(DBI)
library(fixest)
con_fe <- dbConnect(duckdb::duckdb(config=list(temp_directory="")))
set.seed(8721)
d_fe <- data.frame(id=1:720,x=rnorm(720),z=rnorm(720),a=sample(1:3,720,TRUE),o=runif(720,-.2,.2))
for(j in 1:6) d_fe[[paste0("f",j)]] <- sample(1:(3+j),720,TRUE)
d_fe$cluster <- sample(1:18,720,TRUE)
d_fe$y <- rbinom(720,1,plogis(-.7+.3*d_fe$x-.2*d_fe$z+.12*d_fe$f1))
source(system.file('tinytest/helper_glm_reference.R',package='dbreg'),local=TRUE)
reference_fe <- dbglm_reference_feglm
check_fe <- function(actual,reference,atol,rtol) {
 expect_identical(names(actual),names(reference))
 expect_true(all(is.finite(actual)))
 expect_true(all(abs(actual-reference)<=atol+rtol*abs(reference)),
             info=paste(case_fe,"max error",max(abs(actual-reference))))
}
for(link in c("logit","cloglog")) for(nfe in c(1,2,6)) for(kind in c("iid","hc1","cluster")) {
 case_fe <- paste(link,nfe,kind)
 fml <- as.formula(paste("y~x+z|",paste0("f",1:nfe,collapse="+")))
 vc <- if(kind=="cluster") ~f1 else kind
 vr <- if(kind=="hc1") "hetero" else vc
 m <- dbglm(fml,conn=con_fe,data=d_fe,family=binomial(link),weights="a",offset="o",vcov=vc)
 r <- reference_fe(fml,d_fe,family=binomial(link),weights=~a,offset=~o,vcov=vr)
 check_fe(coef(m),coef(r),1e-7,1e-6)
 check_fe(sqrt(diag(vcov(m))),sqrt(diag(vcov(r))),1e-8,1e-6)
 check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
 expect_equal(as.numeric(logLik(m)),as.numeric(logLik(r)),tolerance=1e-8)
 expect_equal(m$n_params_model,r$nparams)
 expect_equal(m$n_params_vcov,attr(vcov(r,attr=TRUE),"df.K"))
 expect_equal(unname(m$n_fe_levels),unname(r$fixef_sizes))
 expect_true(m$score_max<=m$control$score_tol)
 expect_true(m$ap_residual<=m$control$fe_tol)
 expect_true(all(diff(m$deviance_history)<1e-9))
 expect_equal(length(dbListTables(con_fe)),0L)
}
# Complete subgroup slopes, including numeric factor calls and continuous reuse.
for(link in c("cloglog","logit")) {
 case_fe <- paste(link,"full subgroup")
 fml <- y~x:factor(f1)+z+f1|f1+f2
 m <- dbglm(fml,conn=con_fe,data=d_fe,family=binomial(link),vcov=~cluster)
 r <- reference_fe(fml,d_fe,family=binomial(link),vcov=~cluster)
 # The continuous f1 column is absorbed; all four wage slopes must remain.
 expect_equal(length(grep("x:factor(f1)",names(coef(m)),fixed=TRUE)),4L)
 check_fe(coef(m),coef(r),1e-7,1e-6)
 check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
 expect_true("f1" %in% m$collin.var)
}
# Poor starts are damped inside the additive model space.
for(link in c("cloglog","logit")) {
 good <- dbglm(y~x+z|f1+f2,conn=con_fe,data=d_fe,family=binomial(link))
 bad <- dbglm(y~x+z|f1+f2,conn=con_fe,data=d_fe,family=binomial(link),start=c(8,-4))
 expect_true(bad$step_halvings>0)
 expect_equal(coef(bad),coef(good),tolerance=1e-7)
 expect_equal(vcov(bad),vcov(good),tolerance=1e-7)
 expect_true(all(diff(bad$deviance_history)<1e-9))
}
# Rare mixed-outcome FE cells can make a full Fisher step nearly deterministic.
# Backtracking must retain a valid finite fit rather than collapse information.
set.seed(609)
rare <- data.frame(x=rnorm(802),f=sample(1:4,802,TRUE),g=sample(1:3,802,TRUE))
rare$y <- rbinom(802,1,binomial("cloglog")$linkinv(-2+.6*rare$x))
rare$f[801:802] <- 5;rare$x[801:802] <- c(1.2,-1.2);rare$y[801:802] <- c(0,1)
case_fe <- "rare mixed cell"
m <- dbglm(family = stats::binomial("cloglog"), y~x|f+g,conn=con_fe,data=rare)
r <- reference_fe(y~x|f+g,rare,family=binomial("cloglog"))
check_fe(coef(m),coef(r),1e-7,1e-6)
check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
expect_true(m$score_max<=m$control$score_tol)
# No slopes; exact likelihood and inference dimensions still defined.
case_fe <- "FE-only"
m <- dbglm(family = stats::binomial("cloglog"), y~1|f1+f2,conn=con_fe,data=d_fe)
r <- reference_fe(y~1|f1+f2,d_fe,family=binomial("cloglog"))
expect_equal(length(coef(m)),0L)
expect_equal(dim(vcov(m)),c(0L,0L))
expect_equal(as.numeric(logLik(m)),as.numeric(logLik(r)),tolerance=1e-8)
expect_equal(m$n_params_model,r$nparams)
expect_equal(nrow(tidy(m)),0L)
expect_error(predict(m,d_fe),"absorbed")
# Composite FE and multiple nested effects; model df independent of covariance.
d_fe$fsmall <- d_fe$f1 %% 2
d_fe$external <- d_fe$fsmall
d_fe$allcluster <- interaction(d_fe$f1,d_fe$f2,drop=TRUE)
for(cl in c("cluster","f1","fsmall","external","allcluster")) for(policy in c("none","full","nonnested")) {
 case_fe <- paste("corrections",cl,policy)
 vc <- as.formula(paste0("~",cl))
 fml <- y~x|f1+fsmall+f2
 m <- dbglm(family = stats::binomial("cloglog"), fml,conn=con_fe,data=d_fe,vcov=vc,ssc=list(K.fixef=policy))
 r <- reference_fe(fml,d_fe,family=binomial("cloglog"),vcov=vc,ssc=fixest::ssc(K.fixef=policy))
 check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
 expect_equal(m$n_params_model,r$nparams)
 expect_equal(m$n_params_vcov,attr(vcov(r,attr=TRUE),"df.K"))
}
case_fe <- "joint FE"
m <- dbglm(family = stats::binomial("cloglog"), y~x|f1^f2+f3,conn=con_fe,data=d_fe,vcov=~f1)
r <- reference_fe(y~x|f1^f2+f3,d_fe,family=binomial("cloglog"),vcov=~f1)
check_fe(coef(m),coef(r),1e-7,1e-6)
check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
case_fe <- "all nested multiple effects"
m <- dbglm(family = stats::binomial("cloglog"), y~x|f1+fsmall,conn=con_fe,data=d_fe,vcov=~external)
r <- reference_fe(y~x|f1+fsmall,d_fe,family=binomial("cloglog"),vcov=~external)
check_fe(as.vector(vcov(m)),as.vector(vcov(r)),1e-8,1e-5)
expect_equal(m$n_params_vcov,attr(vcov(r,attr=TRUE),"df.K"))
case_fe <- "literal duplicate FE"
m <- dbglm(family = stats::binomial("cloglog"), y~x|f1+f1,conn=con_fe,data=d_fe)
r <- reference_fe(y~x|f1+f1,d_fe,family=binomial("cloglog"))
check_fe(coef(m),coef(r),1e-7,1e-6)
expect_equal(m$n_params_model,r$nparams)
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f1+f2,conn=con_fe,data=d_fe,control=dbglm_control(fe_max_iter=1,fe_tol=1e-30)),"projection")
expect_equal(length(dbListTables(con_fe)),0L)
dbDisconnect(con_fe,shutdown=TRUE)
