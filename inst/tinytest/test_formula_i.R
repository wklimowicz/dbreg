source(system.file("tinytest/helper_glm_reference.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_glm_validation.R",package="dbreg"),local=TRUE)
source(system.file("tinytest/helper_formula_i.R",package="dbreg"),local=TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  dbreg_test_formula_i(conn)
  restore_i = dbreg_reference_i()
  on.exit(restore_i(),add=TRUE)
  i = fixest::i
  set.seed(713)
  d = data.frame(f=rep(letters[1:3],80),j=rep(1:3,80),x=rnorm(240),z=rnorm(240),g=rep(1:8,30),id=1:240)
  d[["y"]] = .3*d[["x"]]+as.numeric(factor(d[["f"]]))+rnorm(240)
  for (fml in list(y~0+i(f),y~0+i(f)+i(g),y~x+i(f,x),y~i(f,ref=TRUE),
    y~i(f,ref="@^[ab]"),y~i(f,keep="B"),y~i(f)+i(f),y~i(f)+factor(f)+z,
    y~i(f,x,keep=c("a","c")),y~j+i(j))) {
    fit = dbreg::feols(fml,d)
    ref = do.call(fixest::feols,list(fml=fml,data=d,notes=FALSE))
    expect_equal(coef(fit),coef(ref),tolerance=1e-7)
    expect_equal(predict(fit,d),predict(ref,d),tolerance=1e-7)
    close(fit)
  }
  # Installed fixest 0.14.0 overflows its stack on explicit ref=NULL. The
  # unreferenced indicator basis is nevertheless exactly the lm dummy space.
  fit = dbreg::feols(y~i(f,ref=NULL),d,lean=TRUE)
  ref = lm(y~f,d)
  expect_equal(predict(fit,d),as.numeric(predict(ref,d)),tolerance=1e-8)
  expect_identical(names(coef(fit)),c("(Intercept)","f::a","f::b"))
  saved = tempfile(fileext=".rds")
  on.exit(unlink(saved),add=TRUE)
  saveRDS(fit,saved)
  close(fit)
  restored = readRDS(saved)
  expect_equal(predict(restored,d),as.numeric(predict(ref,d)),tolerance=1e-8)
  expect_true(as.numeric(object.size(restored))<150000)
  for (args in list(list(d[["f"]]),list(d[["f"]],d[["x"]]),list(d[["f"]],ref="b"),
    list(d[["j"]],ref=1),list(d[["f"]],keep="@^[ac]"))) {
    expected = do.call("i",args)
    actual = do.call(dbreg::i,args)
    expect_equal(actual,expected,tolerance=1e-12)
  }
  expect_error(dbreg::feols(y~i(f,f),d),"second i.*numeric.*logical")
  expect_error(dbreg::feols(y~i(f,bin=2),d),"not supported yet: bin")
  expect_error(dbreg::feols(y~i(f,ref=z),d),"Could not evaluate.*ref")
  expect_error(dbreg::feols(y~i(f):x,d),"additive")
  expect_error(dbreg::feols(y~i(f):x,conn=conn,table="absent"),"additive")
  expect_error(dbreg::feols(y~i(f)*x,conn=conn,table="absent"),"additive")
  expect_error(dbreg::feols(y~x|i(f),d),"fixed-effect block")
  expect_error(dbreg::feols(y~x|g,conn=conn,table="absent",strategy="mundlak",lean=TRUE),"Mundlak.*native dbreg")
  ambiguous = d
  ambiguous[["f"]] = rep(c("a","a:x","b"),80)
  expect_error(dbreg::feols(y~i(f)+i(f,x),ambiguous),"duplicate public coefficient names")
  unsafe = d
  unsafe[["group name"]] = unsafe[["f"]]
  named = dbreg::feols(y~i(`group name`),unsafe)
  named_ref = fixest::feols(y~i(`group name`),unsafe,notes=FALSE)
  expect_equal(coef(named),coef(named_ref),tolerance=1e-8)
  expect_equal(predict(named,unsafe),predict(named_ref,unsafe),tolerance=1e-8)
  output = predict(named,unsafe,conn=conn,row_id="id")
  value = as.data.frame(output)
  expect_equal(value[[output[["prediction_column"]]]][match(unsafe[["id"]],value[["id"]])],predict(named_ref,unsafe),tolerance=1e-8)
  close(output)
  close(named)
  fit = dbreg::feols(y~i(f,x),d)
  unseen = d
  unseen[["f"]][1] = "new"
  expect_error(predict(fit,unseen),"unseen i")
  expect_error(predict(fit,unseen,conn=conn,row_id="id"),"unseen")
  close(fit)
  before = DBI::dbListTables(conn)
  old = options(dbreg.i_max_columns=2)
  on.exit(options(old),add=TRUE)
  expect_error(dbreg::feols(y~i(f,x),d,conn=conn),"i.*columns exceed")
  expect_identical(DBI::dbListTables(conn),before)
  options(old)
  old = options(dbreg.i_max_columns=0)
  expect_error(dbreg::feols(y~i(f),d),"positive integer")
  plain = dbreg::feglm(y~x,d)
  expect_true(all(is.finite(coef(plain))))
  close(plain)
  options(old)
  empty = dbreg::feols(y~i(f,ref=c("a","b","c")),d)
  expect_true(grepl("Database OLS",capture.output(print(empty))[[1L]],fixed=TRUE))
  expect_equal(unname(coef(empty)),mean(d[["y"]]),tolerance=1e-8)
  expect_equal(predict(empty,d),rep(mean(d[["y"]]),nrow(d)),tolerance=1e-8)
  close(empty)
  d[["b"]] = rbinom(nrow(d),1,pnorm(.2*d[["x"]]))
  d[["p"]] = rpois(nrow(d),exp(.2*d[["x"]]))
  for (family in list(gaussian(),poisson(),binomial("logit"),binomial("probit"))) {
    fml = if(family[["family"]]=="gaussian") y~x+i(f,z,ref="a")|g else
      if(family[["family"]]=="poisson") p~x+i(f,z,ref="a")|g else b~x+i(f,z,ref="a")|g
    fit = dbreg::feglm(fml,d,family=family,notes=FALSE)
    ref = dbglm_reference_feglm(fml,d,family=family)
    expect_equal(coef(fit),coef(ref),tolerance=1e-6)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-6)
    expect_equal(predict(fit,d),predict(ref,d),tolerance=1e-6)
    close(fit)
  }
})
