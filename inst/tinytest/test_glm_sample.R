library(dbreg)
library(fixest)
library(DBI)
cs <- dbConnect(duckdb::duckdb(config=list(temp_directory="")))
set.seed(584)
ds <- data.frame(id=1:480,x=rnorm(480),f=rep(1:4,120),g=rep(1:6,each=80),a=1,o=0,cl=rep(1:12,each=40))
ds$y <- rbinom(480,1,plogis(-.5+.2*ds$x))
extra <- ds[1:4,];extra$id <- 481:484;extra$f <- c(5,5,6,6);extra$g <- c(1,7,7,8);extra$y <- c(0,1,0,1)
ds <- rbind(ds,extra)
ds$x[1] <- NA;ds$o[2] <- Inf;ds$a[3] <- 0;ds$cl[4] <- NA
for(mode in c("perfect_fit","singletons","infinite_coef")) {
 m <- dbglm(family = stats::binomial("cloglog"), y~x|f+g,conn=cs,data=ds,weights="a",offset="o",vcov=~cl,
             row_id="id",sample_table="retained",fixef_rm=mode)
 eligible <- ds[-(1:4),]
 r <- fixest::feglm(y~x|f+g,eligible,family=binomial("cloglog"),weights=~a,offset=~o,vcov=~cl,
             fixef.rm=mode,glm.tol=3e-13,fixef.tol=1e-10,notes=FALSE)
 keys <- sort(dbReadTable(cs,"retained")$id)
 expect_equal(keys,sort(eligible$id[fixest::obs(r)]))
 expect_equal(nobs(m),nobs(r))
 expect_true(m$pruning_rounds>1)
 expect_equal(sum(m$exclusions),4)
 expect_equal(m$nobs_dropped,8)
 expect_equal(coef(m),coef(r),tolerance=1e-7)
 expect_equal(dbListTables(cs),"retained")
 expect_error(dbglm(family = stats::binomial("cloglog"), y~x,conn=cs,data=ds,sample_table="retained"),"exists")
 dbRemoveTable(cs,"retained")
}
# All modes coincide when there are no perfect groups; none preserves all keys.
clean <- ds[5:480,]
for(mode in c("none","singletons","infinite_coef","perfect_fit")) {
 m <- dbglm(family = stats::binomial("cloglog"), y~x|f+g,conn=cs,data=clean,fixef_rm=mode,row_id=c("id","f"),sample_table="retained")
 expect_equal(sort(dbReadTable(cs,"retained")$id),sort(clean$id))
 expect_equal(m$nobs_pruned,0)
 dbRemoveTable(cs,"retained")
}
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f,conn=cs,data=transform(clean,id=1),row_id="id"),"unique")
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f,conn=cs,data=transform(clean,id=NA),row_id="id"),"nonmissing")
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f,conn=cs,data=clean,control=dbglm_control(glm_max_iter=1),sample_table="failed"),"converge")
expect_false(dbExistsTable(cs,"failed"))
expect_equal(length(dbListTables(cs)),0L)
# Missing outcome groups cannot return ordinary successful inference under none.
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f+g,conn=cs,data=ds,fixef_rm="none"),"boundary|converge|singular|halving|separation")
expect_equal(length(dbListTables(cs)),0L)
# Type metadata for factor(column) is distinct from continuous column use.
clean$code <- rep(c(1,2,10,20),length.out=nrow(clean))
m <- dbglm(y~x+factor(code)+code,conn=cs,data=clean,family=binomial("logit"),ssc=list(K.adj=FALSE))
r <- glm(y~x+factor(code),clean,family=binomial("logit"),control=glm.control(epsilon=1e-14))
expect_true("code" %in% m$collin.var)
expect_equal(coef(m),coef(r),tolerance=1e-7)
expect_equal(predict(m,clean[1,],type="response"),unname(predict(r,clean[1,],type="response")),tolerance=1e-7)
expect_equal(as.numeric(m$model_spec$factor_levels[[1]]),c(1,2,10,20))
expect_error(predict(m,transform(clean[1,],code=99)),"unseen")
expect_equal(coef(dbglm(family = stats::binomial("cloglog"), y~x,conn=cs,data=transform(clean,a=.5),weights="a")),
             coef(dbglm(family = stats::binomial("cloglog"), y~x,conn=cs,data=clean)),tolerance=1e-7)
expect_equal(length(dbListTables(cs)),0L)
# Wrapper calls stay compact even if do.call supplies literal rows/connections.
big <- data.frame(x=seq(-2,2,length.out=2000),y=rep(c(0,1),1000),
                  unused=paste0(seq_len(2000),strrep("x",200)))
ff <- local({ captured <- big; y~x })
compact <- do.call(dbglm,list(fml=ff,conn=cs,data=big,family=binomial('logit')))
expect_true(is.symbol(compact$call$data))
expect_true(is.symbol(compact$call$conn))
expect_identical(environment(compact$fml),baseenv())
expect_identical(environment(compact$model_spec$rhs),baseenv())
expect_true(length(serialize(compact,NULL))<100000)
expect_equal(coef(compact),coef(dbglm("y~x",conn=cs,data=big,family=binomial("logit"))))
expect_equal(length(dbListTables(cs)),0L)
# Invalid outcomes cannot disappear unnoticed through perfect-group pruning.
invalid <- clean
invalid$y[invalid$f==max(invalid$f)] <- 2
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f,conn=cs,data=invalid),"zero or one")
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f|g,conn=cs,data=clean),"at most")
expect_error(dbglm(family = stats::binomial("cloglog"), y~x|f|x~g,conn=cs,data=clean),"outcome|at most")
expect_error(dbglm(family = stats::binomial("cloglog"), y~I(x^2)|f,conn=cs,data=transform(clean,x=1e308)),"nonfinite transformed")
expect_equal(length(dbListTables(cs)),0L)
dbDisconnect(cs,shutdown=TRUE)
