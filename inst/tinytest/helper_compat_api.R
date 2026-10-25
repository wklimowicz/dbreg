dbreg_test_compat_api = function(conn) {
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  models = outputs = list()
  on.exit({ for (output in outputs) close(output); for (model in models) close(model) }, add = TRUE)
  set.seed(582)
  n = 360L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:12, each = 30),
    h = sample(1:6,n,TRUE), f = rep(c("a","b","c"),120), o = rnorm(n,sd=.1),
    w = runif(n,.1,1.3), keep = rep(c(TRUE,TRUE,FALSE,TRUE),90), u = 1)
  d[["cl"]] = d[["h"]]
  d[["y"]] = .4*d[["x"]]-.2*d[["z"]]+sin(d[["g"]])+.1*cos(d[["h"]])+d[["o"]]+rnorm(n,sd=.4)
  d[["b"]] = rbinom(n,1,1-exp(-exp(.2*d[["x"]]+.1*sin(d[["g"]])+d[["o"]])))
  d[["x"]][2] = NA_real_
  d[["o"]][3] = NA_real_
  d[["w"]][5] = 0
  d[["cl"]][6] = NA_integer_
  d[["keep"]][8] = NA
  d[["u"]][17:18] = c(NA_real_,NaN)
  d = rbind(d,d[9,])
  d[["id"]][n+1L] = n+1L
  d[["g"]][n+1L] = 99L
  d[["__dbglm_single"]] = TRUE
  d[["__dbglm_infinite"]] = "unused"
  table = dbreg:::prediction_upload_frame(conn,d,resources)
  # Prediction uploads intentionally treat R NaNs as missing. Preserve a true
  # numeric SQL NaN here to exercise predicate semantics independently.
  raw = table
  table = basename(tempfile("dbreg_subset_nan_"))
  resources[["tables"]] = c(resources[["tables"]],table)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn,x))
  query = paste("SELECT",paste(c(q(setdiff(names(d),"u")),
    "CASE WHEN id=18 THEN CAST('NaN' AS DOUBLE) ELSE u END AS u"),collapse=","),"FROM",raw)
  DBI::dbExecute(conn,dbreg:::sql_sample_table(conn,table,query,names(d)))
  baseline = dbglm_validation_tables(conn)
  cases = list(
    list(estimator="feols",fml=y~x+z|g+h,vcov=~g+cl),
    list(estimator="feols",fml=y~x|g[z],vcov=~g+cl),
    list(estimator="feols",fml=y~x+factor(f)|g,vcov="hetero"),
    list(estimator="feglm",fml=b~x+z|g+h,family=binomial("cloglog"),vcov=~cl),
    list(estimator="feglm",fml=y~x+z|g,family="gaussian",vcov="iid"),
    list(estimator="feols",fml=y~x+z,vcov=~cl,no_subset=TRUE))
  for (case in cases) {
    message("[Compatibility API] ",case[["estimator"]]," ",paste(deparse(case[["fml"]]),collapse=""))
    args = list(fml=case[["fml"]],conn=conn,table=table,weights=~w,offset=~o,
      subset=~keep & u>0,vcov=case[["vcov"]],row_id="id",notes=FALSE)
    refargs = list(fml=case[["fml"]],data=d,weights=~w,offset=~o,subset=~keep & u>0,
      vcov=case[["vcov"]],notes=FALSE)
    if (isTRUE(case[["no_subset"]])) args[["subset"]] = refargs[["subset"]] = NULL
    if (case[["estimator"]] == "feglm") {
      args[["family"]] = refargs[["family"]] = case[["family"]]
      refargs[["glm.tol"]] = 1e-12
      refargs[["fixef.tol"]] = 1e-10
    }
    fit = do.call(if(case[["estimator"]]=="feols") dbreg::feols else dbreg::feglm,args)
    models[[length(models)+1L]] = fit
    ref = do.call(if(case[["estimator"]]=="feols") fixest::feols else dbglm_reference_feglm,refargs)
    expect_equal(coef(fit),coef(ref),tolerance=1e-6)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-6)
    expect_equal(as.numeric(dbreg::coeftable(fit)),as.numeric(as.matrix(fixest::coeftable(ref))),tolerance=1e-5)
    expect_equal(as.numeric(fixest::se(fit)),as.numeric(fixest::se(ref)),tolerance=1e-6)
    expect_equal(as.numeric(fixest::pvalue(fit)),as.numeric(fixest::pvalue(ref)),tolerance=1e-5)
    expect_identical(names(fixest::se(fit)),names(coef(ref)))
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-6)
    expect_equal(AIC(fit),AIC(ref),tolerance=1e-6)
    expect_equal(BIC(fit),BIC(ref),tolerance=1e-6)
    expect_equal(df.residual(fit),df.residual(ref))
    expect_equal(dbreg::degrees_freedom(fit,"k"),fixest::degrees_freedom(ref,"k"))
    expect_equal(dbreg::degrees_freedom(fit,"t"),fixest::degrees_freedom(ref,"t"))
    ids = fixest::obs(ref)
    output = dbreg::obs(fit)
    values = as.data.frame(output)
    expect_equal(sort(values[["id"]]),ids)
    close(output)
    output = predict(fit)
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(ids,values[["id"]])],predict(ref),tolerance=1e-6)
    close(output)
    expect_equal(predict(fit,d),predict(ref,d),tolerance=1e-6)
    expect_equal(coef(fit,keep="x",order="x"),coef(ref,keep="x",order="x"),tolerance=1e-6)
    close(fit)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-6)
    expect_identical(dbglm_validation_tables(conn),baseline)
  }
  # Predicate-only NaNs must follow R comparisons, and %in% distinguishes NA
  # from NaN even though is.na identifies both.
  schema = DBI::dbGetQuery(conn,paste("SELECT * FROM",table,"LIMIT 0"))
  for (predicate in list(~u>0,~u %in% c(NA,1),~u %in% c(NaN,1),~f %in% c("a","c"),~!(f=="b"),~keep==TRUE)) {
    where = dbreg:::sql_subset_predicate(conn,predicate,schema)
    actual = dbreg:::backend_get_query(conn,paste("SELECT id FROM",table,"WHERE",where))[["id"]]
    expected = eval(predicate[[2L]],d,baseenv())
    expected[is.na(expected)] = FALSE
    expect_equal(sort(actual),d[["id"]][expected])
  }
  fit = dbreg(y~x+z|g+h,conn=conn,table=table,weights="w",offset="o",subset=~keep & u>0,
    fixef_rm="singletons",retain="inference",row_id="id")
  models[[length(models)+1L]] = fit
  ref = fixest::feols(y~x+z|g+h,d,weights=~w,offset=~o,subset=~keep & u>0,notes=FALSE)
  expect_identical(fit[["fitted_state"]][["prediction"]],NULL)
  expect_equal(as.numeric(vcov(fit,vcov=~g+h,ssc=dbreg::ssc())),as.numeric(vcov(ref,vcov=~g+h)),tolerance=1e-6)
  output = predict(fit)
  values = as.data.frame(output)
  expect_equal(values[[output[["prediction_column"]]]][match(fixest::obs(ref),values[["id"]])],predict(ref),tolerance=1e-6)
  close(output)
  close(fit)
  expect_identical(dbglm_validation_tables(conn),baseline)
}
