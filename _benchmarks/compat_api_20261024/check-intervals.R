pkgload::load_all()
source('inst/tinytest/helper_glm_reference.R')
set.seed(585)
d=data.frame(x=rnorm(240),g=rep(1:12,20),h=rep(1:8,30),w=runif(240,.2,1.3))
d[['y']]=d[['x']]+sin(d[['g']])+rnorm(240)
d[['b']]=rbinom(240,1,plogis(.3*d[['x']]))
for(engine in c('feols','feglm')) {
 if(engine=='feols') {
  fit=dbreg::feols(y~x|g+h,d,weights=~w,vcov=~g+h,notes=FALSE)
  ref=fixest::feols(y~x|g+h,d,weights=~w,vcov=~g+h,notes=FALSE)
 } else {
  fit=dbreg::feglm(b~x|g,d,weights=~w,vcov=~h,family=binomial('cloglog'),notes=FALSE)
  ref=suppressWarnings(dbglm_reference_feglm(b~x|g,d,weights=~w,vcov=~h,family=binomial('cloglog')))
 }
 for(policy in list(NULL,dbreg::ssc(t.df=7))) {
  actual=confint(fit,ssc=policy)
  expected=as.matrix(confint(ref,ssc=policy))
  error=max(abs(actual-expected))
  stopifnot(error<1e-6)
  cat(engine,if(is.null(policy)) 'default' else 'numeric t.df','max interval error:',error,'\n')
 }
 close(fit)
}
