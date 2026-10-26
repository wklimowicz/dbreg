# Every assertion compares observable values/coding with fixest or analytical
# dummy regressions; both engines run the same formula/prediction fixture.
dbreg_test_formula_i = function(conn) {
  predict_fit = function(object, ...) {
    if (inherits(object,"dbglm")) predict(object,...,type="response") else predict(object,...)
  }
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources),add=TRUE)
  models = outputs = list()
  on.exit({ for(output in outputs) close(output); for(model in models) close(model) },add=TRUE)
  i = fixest::i
  set.seed(712)
  n = 240L
  d = data.frame(id=seq_len(n),f=rep(c("alpha","Beta","dbglm_factor_i1__quote'\\slash"),80),
    j=rep(c(1L,2L,10L),80),flag=rep(c(FALSE,TRUE),120),x=rnorm(n),z=rnorm(n),
    g=rep(1:8,each=30),h=sample(1:6,n,TRUE),w=runif(n,.2,1.4),o=rnorm(n,sd=.1))
  d[["y"]] = .3*d[["x"]]+.2*d[["z"]]+.2*(d[["j"]]==2)+sin(d[["g"]])+.2*rnorm(n)+d[["o"]]
  d[["b"]] = rbinom(n,1,1-exp(-exp(.2*d[["x"]]+.1*(d[["j"]]==2)+d[["o"]])))
  d[["y"]][1] = NA_real_
  d[["x"]][3] = NA_real_
  d[["w"]][4] = 0
  table = dbreg:::prediction_upload_frame(conn,d,resources)
  baseline = dbglm_validation_tables(conn)
  cases = list(
    list(fml=y~x+i(f)),list(fml=y~i(f,x)),list(fml=y~x+i(f,x,ref="Beta")|g),
    list(fml=y~i(f,keep=c("alpha","Beta"))),list(fml=y~i(j,2)),
    list(fml=y~i(flag)),list(fml=y~x+i(f,flag)|g),
    list(fml=b~x+i(f)|g,family=binomial("cloglog")),
    list(fml=y~x+i(f,z,ref="alpha")|g+h),
    list(fml=y~i(f,x)|g[z]))
  for (case in cases) {
    fml = case[["fml"]]
    glm = !is.null(case[["family"]])
    message("[i() formula] ",paste(deparse(fml),collapse=""))
    args = list(fml=fml,conn=conn,table=table,weights=~w,offset=~o,row_id="id",retain="all",ssc=fixest::ssc(),fixef_rm="none",
      vcov=if(glm) ~g else ~g+h)
    refargs = list(fml=fml,data=d,weights=~w,offset=~o,notes=FALSE,vcov=args[["vcov"]])
    if(glm) args[["family"]] = refargs[["family"]] = case[["family"]]
    requested_vcov = args[["vcov"]]
    args[["vcov"]] = if(inherits(requested_vcov,"formula")) stats::reformulate(tail(all.vars(requested_vcov),1L)) else requested_vcov
    fit = do.call(if(glm) dbreg::dbglm else dbreg::dbreg,args)
    fit = summary(fit,vcov=requested_vcov)
    models[[length(models)+1L]] = fit
    ref = do.call(if(glm) dbglm_reference_feglm else fixest::feols,refargs)
    expect_equal(coef(fit),coef(ref),tolerance=1e-6)
    expect_equal(as.numeric(vcov(fit)),as.numeric(vcov(ref)),tolerance=1e-6)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-6)
    expect_equal(predict_fit(fit,d),predict(ref,d),tolerance=1e-6)
    expect_identical(names(coef(fit)),names(coef(ref)))
    if (identical(case,cases[[3L]])) {
      mixed = d[10:12,]
      mixed[["f"]][1:2] = "unknown"
      mixed[["g"]][1] = 999L
      mixed[["x"]][2] = NA_real_
      expected = c(NA_real_,NA_real_,unname(predict(ref,mixed[3,,drop=FALSE])))
      expect_equal(predict_fit(fit,mixed),expected,tolerance=1e-6)
      output = predict_fit(fit,mixed,conn=conn,row_id="id")
      outputs[[length(outputs)+1L]] = output
      value = as.data.frame(output)
      expect_equal(value[[output[["prediction_column"]]]][match(mixed[["id"]],value[["id"]])],expected,tolerance=1e-6)
      close(output)
    }
    for (frame in c(FALSE,TRUE)) {
      output = if(frame) predict_fit(fit,d,conn=conn,row_id="id") else predict_fit(fit,conn=conn,table=table,row_id="id")
      outputs[[length(outputs)+1L]] = output
      value = as.data.frame(output)
      value = value[[output[["prediction_column"]]]][match(d[["id"]],value[["id"]])]
      expect_equal(value,predict(ref,d),tolerance=1e-6)
      close(output)
    }
    output = predict_fit(fit)
    outputs[[length(outputs)+1L]] = output
    value = as.data.frame(output)
    ids = fixest::obs(ref)
    expect_equal(value[[output[["prediction_column"]]]][match(ids,value[["id"]])],predict(ref),tolerance=1e-6)
    close(output)
    expect_equal(as.numeric(vcov(fit,vcov="hetero")),as.numeric(vcov(ref,vcov="hetero")),tolerance=1e-6)
    close(fit)
    expect_identical(dbglm_validation_tables(conn),baseline)
  }
  # Driver conversion and display rounding must not merge fitted categories.
  for (kind in c("integer","numeric")) {
    keys = if(kind=="integer") c("9007199254740993","9007199254740994","9007199254740995") else 1+c(0,1,2)*.Machine[["double.eps"]]
    key_data = data.frame(id=1:90,x=rnorm(90),key=rep(keys,30),index=factor(rep(1:3,30)))
    key_data[["y"]] = key_data[["x"]]+as.numeric(key_data[["index"]])+rnorm(90,sd=.1)
    raw = dbreg:::prediction_upload_frame(conn,key_data,resources)
    typed = raw
    if(kind=="integer") {
      typed = basename(tempfile("dbreg_i_exact_"))
      resources[["tables"]] = c(resources[["tables"]],typed)
      query = paste("SELECT id,x,y,CAST(key AS BIGINT) AS key FROM",raw)
      DBI::dbExecute(conn,dbreg:::sql_sample_table(conn,typed,query,c("id","x","y","key")))
    }
    fit = dbreg::dbreg(y~x+i(key),conn=conn,table=typed,row_id="id",retain="all")
    models[[length(models)+1L]] = fit
    ref = lm(y~x+index,key_data)
    truth = as.numeric(predict(ref,key_data))
    expect_equal(predict_fit(fit,key_data),truth,tolerance=1e-8)
    expect_equal(unname(coef(fit)),unname(coef(ref)),tolerance=1e-8)
    expect_equal(length(fit[["i_specs"]][[1]][["keys"]]),3L)
    for (frame in c(FALSE,TRUE)) {
      output = if(frame) predict_fit(fit,key_data,conn=conn,row_id="id") else predict_fit(fit,conn=conn,table=typed,row_id="id")
      outputs[[length(outputs)+1L]] = output
      value = as.data.frame(output)
      expect_equal(value[[output[["prediction_column"]]]][match(key_data[["id"]],value[["id"]])],truth,tolerance=1e-8)
      close(output)
    }
    if(kind=="integer") {
      lossy = key_data
      lossy[["key"]] = as.numeric(lossy[["key"]])
      expect_error(predict_fit(fit,lossy),"must be exact")
      expect_error(predict_fit(fit,lossy,conn=conn,row_id="id"),"must be exact")
      malformed = key_data
      malformed[["key"]][1] = "not-an-integer"
      expect_error(predict_fit(fit,malformed),"unseen i")
      expect_error(predict_fit(fit,malformed,conn=conn,row_id="id"),"unseen")
      incomplete = key_data[1,,drop=FALSE]
      incomplete[["key"]] = as.numeric(incomplete[["key"]])
      incomplete[["x"]] = NA_real_
      expect_true(is.na(predict_fit(fit,incomplete)))
      output = predict_fit(fit,incomplete,conn=conn,row_id="id")
      outputs[[length(outputs)+1L]] = output
      value = as.data.frame(output)
      expect_true(is.na(value[[output[["prediction_column"]]]]))
      close(output)
    } else {
      missing = key_data
      missing[["key"]][1:2] = c(NA_real_,Inf)
      expected = truth
      expected[1:2] = NA_real_
      expect_equal(predict_fit(fit,missing),expected,tolerance=1e-8)
      output = predict_fit(fit,missing,conn=conn,row_id="id")
      outputs[[length(outputs)+1L]] = output
      value = as.data.frame(output)
      expect_equal(value[[output[["prediction_column"]]]][match(missing[["id"]],value[["id"]])],expected,tolerance=1e-8)
      close(output)
    }
    close(fit)
  }
}
