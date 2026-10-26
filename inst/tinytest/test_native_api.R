source(system.file("tinytest/helper_glm_reference.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_glm_validation.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_native_api.R",package="dbreg"),local=TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  dbreg_test_native_api(conn)
  set.seed(583)
  d = data.frame(x=rnorm(160),z=rnorm(160),g=rep(1:8,20),o=rnorm(160),w=runif(160),keep=rep(c(TRUE,FALSE,TRUE,TRUE),40))
  d[["y"]] = d[["x"]]+d[["o"]]+sin(d[["g"]])+rnorm(160)
  d[["w"]][2] = 0
  d[["y"]][3] = NA_real_
  d[["keep"]][4] = NA
  fit = dbreg::dbreg(y~x+z|g,data=d,weights=~w,offset=~o,subset=~keep,retain="all",ssc=fixest::ssc())
  ref = fixest::feols(y~x+z|g,d,weights=~w,offset=~o,subset=~keep,notes=FALSE)
  expect_identical(stats::case.names(fit),fixest::obs(ref))
  expect_equal(predict(fit),predict(ref),tolerance=1e-7)
  expect_equal(predict(fit,sample="original"),predict(ref,sample="original"),tolerance=1e-7)
  expect_equal(predict(fit,d),predict(ref,d),tolerance=1e-7)
  expect_equal(as.numeric(fixest::coeftable(fit,cluster=~g)),as.numeric(as.matrix(fixest::coeftable(ref,cluster=~g))),tolerance=1e-6)
  for (selector in list(list(keep="!Intercept",order="z"),list(drop="x",order="!z"),list(keep="%x"))) {
    expect_equal(as.numeric(do.call(fixest::coeftable,c(list(fit),selector))),as.numeric(do.call(fixest::coeftable,c(list(ref),selector))),tolerance=1e-7)
  }
  changed = summary(fit,ssc=fixest::ssc(K.fixef="none"))
  changed_ref = summary(ref,ssc=fixest::ssc(K.fixef="none"))
  expect_equal(as.numeric(vcov(changed,cluster=~g)),as.numeric(vcov(changed_ref,cluster=~g)),tolerance=1e-8)
  close(fit)
  expect_error(dbreg::dbreg(y~x,data=d,weights=~log(w)),"one column")
  expect_error(dbreg::dbreg(y~x,data=d,weights=1:5),"column")
  expect_error(dbreg::dbreg(y~x,data=d,subset=c(1,1,2)),"formula")
  d[["b"]] = rbinom(nrow(d),1,plogis(.2*d[["x"]]))
  for (family in list(gaussian(),binomial("cloglog"))) {
    fml = if (family[["family"]]=="gaussian") y~x+z else b~x+z
    fit = dbreg::dbglm(fml,data=d,family=family,weights=~w,offset=~o,subset=~keep,retain="prediction")
    ref = dbglm_reference_feglm(fml,d,family=family,weights=~w,offset=~o,subset=~keep)
    expect_equal(coef(fit),coef(ref),tolerance=1e-6,scale=1)
    expect_equal(predict(fit,d,type="response"),predict(ref,d),tolerance=1e-6)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-6)
    expect_equal(AIC(fit),AIC(ref),tolerance=1e-6)
    expect_identical(fit[["resource_state"]],NULL)
    if (family[["family"]]=="binomial") {
      rows = fixest::obs(ref)
      X = model.matrix(~x+z,d[rows,])
      eta = as.vector(X%*%coef(fit))+d[["o"]][rows]
      mu = family[["linkinv"]](eta)
      score = colSums(X*(d[["w"]][rows]*family[["mu.eta"]](eta)*(d[["b"]][rows]-mu)/(mu*(1-mu))))
      expect_true(max(abs(score))<1e-6)
    }
    close(fit)
  }
  fit = dbreg::dbreg(y~x+z,data=d,offset=~o,retain="prediction")
  saved = tempfile(fileext=".rds")
  on.exit(unlink(saved),add=TRUE)
  saveRDS(fit,saved)
  restored = readRDS(saved)
  expect_equal(predict(restored,d),predict(fit,d))
  expect_identical(logLik(restored),logLik(fit))
  expect_true(as.numeric(object.size(restored))<250000)
  close(fit)
  expect_identical(DBI::dbListTables(conn),character())
})
