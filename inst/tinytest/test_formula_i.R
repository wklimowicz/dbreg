library(fixest)
source(system.file("tinytest/helper_glm_reference.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_glm_validation.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_formula_i.R",package="dbreg"),local=TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  dbreg_test_formula_i(conn)
  i = fixest::i
  set.seed(713)
  d = data.frame(f=rep(letters[1:3],80),j=rep(1:3,80),x=rnorm(240),z=rnorm(240),g=rep(1:8,30),id=1:240)
  d[["y"]] = .3*d[["x"]]+as.numeric(factor(d[["f"]]))+rnorm(240)
  for (fml in list(y~0+i(f),y~0+i(f)+i(g),y~x+i(f,x),y~i(f,ref=TRUE),
    y~i(f,ref="@^[ab]"),y~i(f,keep="B"),y~i(f)+i(f),y~i(f)+factor(f)+z,
    y~i(f,x,keep=c("a","c")),y~j+i(j))) {
    fit = dbreg::dbreg(fml,data=d,retain="all",ssc=fixest::ssc())
    reference_formula = if(identical(fml,y~i(f,ref="@^[ab]"))) y~i(f,ref=c("a","b")) else fml
    ref = do.call(fixest::feols,list(fml=reference_formula,data=d,notes=FALSE))
    expect_equal(coef(fit),coef(ref),tolerance=1e-7)
    expect_equal(predict(fit,d),predict(ref,d),tolerance=1e-7)
    close(fit)
  }
  # Installed fixest 0.14.0 overflows its stack on explicit ref=NULL. The
  # unreferenced indicator basis is nevertheless exactly the lm dummy space.
  fit = dbreg::dbreg(y~i(f,ref=NULL),data=d,retain="prediction")
  ref = lm(y~f,data=d)
  expect_equal(predict(fit,d),as.numeric(predict(ref,d)),tolerance=1e-8)
  expect_identical(names(coef(fit)),c("(Intercept)","f::a","f::b"))
  saved = tempfile(fileext=".rds")
  on.exit(unlink(saved),add=TRUE)
  saveRDS(fit,saved)
  close(fit)
  restored = readRDS(saved)
  expect_equal(predict(restored,d),as.numeric(predict(ref,d)),tolerance=1e-8)
  expect_true(as.numeric(object.size(restored))<150000)
  expect_error(dbreg::dbreg(y~i(f,f),data=d,retain="all",ssc=fixest::ssc()),"second i.*numeric.*logical")
  expect_error(dbreg::dbreg(y~i(f,bin=2),data=d,retain="all",ssc=fixest::ssc()),"not supported yet: bin")
  expect_error(dbreg::dbreg(y~i(f,ref=z),data=d,retain="all",ssc=fixest::ssc()),"Could not evaluate.*ref")
  expect_error(dbreg::dbreg(y~i(f):x,data=d,retain="all",ssc=fixest::ssc()),"interacted")
  expect_error(dbreg::dbreg(y~x|i(f),data=d,retain="all",ssc=fixest::ssc()),"Fixed-effect interactions")
  ambiguous = d
  ambiguous[["f"]] = rep(c("a","a:x","b"),80)
  expect_error(dbreg::dbreg(y~i(f)+i(f,x),data=ambiguous),"duplicate public coefficient names")
  unsafe = d
  unsafe[["group name"]] = unsafe[["f"]]
  named = dbreg::dbreg(y~i(`group name`),data=unsafe,retain="all")
  named_ref = fixest::feols(y~i(`group name`),unsafe,notes=FALSE)
  expect_equal(coef(named),coef(named_ref),tolerance=1e-8)
  expect_equal(predict(named,unsafe),predict(named_ref,unsafe),tolerance=1e-8)
  output = predict(named,unsafe,conn=conn,row_id="id")
  value = as.data.frame(output)
  expect_equal(value[[output[["prediction_column"]]]][match(unsafe[["id"]],value[["id"]])],predict(named_ref,unsafe),tolerance=1e-8)
  close(output)
  close(named)
  fit = dbreg::dbreg(y~i(f,x),data=d,retain="all",ssc=fixest::ssc())
  unseen = d
  unseen[["f"]][1] = "new"
  expect_error(predict(fit,unseen),"unseen i")
  expect_error(predict(fit,unseen,conn=conn,row_id="id"),"unseen")
  close(fit)
  before = DBI::dbListTables(conn)
  old = options(dbreg.i_max_columns=2)
  on.exit(options(old),add=TRUE)
  expect_error(dbreg::dbreg(y~i(f,x),data=d,conn=conn),"i.*columns exceed")
  expect_identical(DBI::dbListTables(conn),before)
  options(old)
  old = options(dbreg.i_max_columns=0)
  expect_error(dbreg::dbreg(y~i(f),data=d,retain="all",ssc=fixest::ssc()),"positive integer")
  plain = dbreg::dbglm(y~x,data=d,retain="all",ssc=fixest::ssc())
  expect_true(all(is.finite(coef(plain))))
  close(plain)
  options(old)
  empty = dbreg::dbreg(y~i(f,ref=c("a","b","c")),data=d,retain="all",ssc=fixest::ssc())
  expect_true(grepl("OLS estimation",capture.output(print(empty))[[1L]],fixed=TRUE))
  expect_equal(unname(coef(empty)),mean(d[["y"]]),tolerance=1e-8)
  expect_equal(predict(empty,d),rep(mean(d[["y"]]),nrow(d)),tolerance=1e-8)
  close(empty)
  d[["b"]] = rbinom(nrow(d),1,pnorm(.2*d[["x"]]))
  d[["p"]] = rpois(nrow(d),exp(.2*d[["x"]]))
  for (family in list(gaussian(),poisson(),binomial("logit"),binomial("probit"))) {
    fml = if(family[["family"]]=="gaussian") y~x+i(f,z,ref="a")|g else
      if(family[["family"]]=="poisson") p~x+i(f,z,ref="a")|g else b~x+i(f,z,ref="a")|g
    fit = dbreg::dbglm(fml,data=d,family=family,retain="all")
    ref = dbglm_reference_feglm(fml,data=d,family=family)
    expect_equal(coef(fit),coef(ref),tolerance=1e-6)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-6)
    expect_equal(predict(fit,d,type="response"),predict(ref,d),tolerance=1e-6)
    close(fit)
  }
})
