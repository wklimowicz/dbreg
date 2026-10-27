# Explicit null diagnostics share the exact point-fit sample and reporting rules.
library(DBI)
local({
  set.seed(9207)
  d = expand.grid(g = 1:6, h = 1:4, replicate = 1:10)
  d[["id"]] = seq_len(nrow(d))
  d[["x"]] = sample(0:3, nrow(d), TRUE)
  d[["w"]] = rep(c(.5,1.5,2.5), length.out = nrow(d))
  d[["o"]] = .05*(d[["x"]] %% 2)
  d[["cl"]] = (d[["id"]] %% 12)+1
  likelihood = function(family, y, mu, w) {
    if (family[["family"]] == "gaussian") {
      n = length(y)
      weighted = -n/2*(log(2*pi)+log(sum(w*(y-mu)^2)/n)+1)+sum(log(w))/2
      return(c(reported = weighted-1, weighted = weighted))
    }
    if (family[["family"]] == "poisson") return(c(reported = sum(w*stats::dpois(y,mu,log=TRUE)),
      weighted = sum(w*stats::dpois(y,mu,log=TRUE))))
    c(reported = sum(round(w)*stats::dbinom(y,1,mu,log=TRUE)), weighted = sum(w*stats::dbinom(y,1,mu,log=TRUE)))
  }
  for (family in list(stats::gaussian(), stats::poisson(), stats::binomial(),
    stats::binomial("probit"), stats::binomial("cloglog"))) {
    eta = -.3+.2*d[["x"]]+d[["o"]]+.1*sin(d[["g"]])
    d[["y"]] = switch(family[["family"]], gaussian = eta+rnorm(nrow(d)),
      poisson = rpois(nrow(d),exp(eta)), binomial = rbinom(nrow(d),1,family[["linkinv"]](eta)))
    fit = dbglm(y ~ x | g+h, data = d, family = family, weights = ~w, offset = ~o,
      vcov = ~cl, fixef_rm = "none", null_model = TRUE)
    record = fit[["null_diagnostics"]]
    reference = suppressWarnings(stats::glm(y ~ 1, data = d, weights = w, offset = o, family = family,
      control = stats::glm.control(epsilon = 1e-12, maxit = 100)))
    expected = likelihood(family, d[["y"]], stats::predict(reference,type="response"), d[["w"]])
    expect_true(record[["fit"]][["available"]])
    expect_equal(record[["nobs"]], nrow(d))
    expect_equal(record[["reported"]][["logLik"]], unname(expected[["reported"]]), tolerance = 1e-7)
    expect_equal(record[["weighted"]][["logLik"]], unname(expected[["weighted"]]), tolerance = 1e-7)
    compressed = dbglm(y ~ x | g+h, data = d, family = family, weights = ~w, offset = ~o,
      vcov = ~cl, fixef_rm = "none", null_model = TRUE, strategy = "compress")
    expect_equal(compressed[["null_diagnostics"]][c("reported","weighted")], record[c("reported","weighted")], tolerance = 1e-7)
    values = gof(fit)
    frame = generics::glance(fit)
    expect_true(is.numeric(values))
    expect_equal(unname(values[["null.logLik"]]), frame[["null.logLik"]])
    expect_equal(unname(values[["pseudo.r.squared"]]), frame[["pseudo.r.squared"]])
    before = capture.output(print(fit))
    close(fit)
    expect_identical(capture.output(print(fit)), before)
    expect_identical(gof(fit), values)
    expect_identical(generics::glance(fit), frame)
    close(compressed)
    # Offset-only null has no constant in a structurally no-intercept model.
    fixed = dbglm(y ~ 0+x, data = d, family = family, weights = ~w, offset = ~o, null_model = TRUE)
    expected_fixed = likelihood(family, d[["y"]], family[["linkinv"]](d[["o"]]), d[["w"]])
    expect_false(fixed[["null_diagnostics"]][["constant"]])
    expect_equal(fixed[["null_diagnostics"]][["reported"]][["logLik"]], unname(expected_fixed[["reported"]]), tolerance = 1e-7)
    # A truly intercept-only full model can reuse its already computed values.
    same = dbglm(y ~ 1, data = d, family = family, weights = ~w, offset = ~o, null_model = TRUE)
    expect_identical(same[["null_diagnostics"]][["fit"]][["method"]], "full_is_null")
    expect_equal(same[["null_diagnostics"]][["reported"]][["ratio"]], 0, tolerance = 0)
    close(fixed); close(same)
  }
  # Binary original outcomes remain binary even though compressed means do not.
  expect_error(dbglm(y ~ x, data = transform(d,y=.25), family = stats::binomial(), null_model = TRUE), "zero or one")
  for (flag in list(NA,1,"yes",logical())) {
    expect_error(dbglm(y ~ x, data = d, null_model = flag), "null_model.*TRUE or FALSE")
  }
  # Reported weights can all round to zero; weighted likelihood remains meaningful.
  d[["y"]] = rep(0:1,length.out=nrow(d)); d[["w"]] = .25
  zero = dbglm(y ~ x, data = d, family = stats::binomial(), weights = ~w, null_model = TRUE)
  expect_identical(zero[["null_diagnostics"]][["reported"]][["ratio_reason"]], "zero_null_likelihood")
  expect_true(zero[["null_diagnostics"]][["weighted"]][["ratio_available"]])
  close(zero)
})

# Default missing/explicit FALSE has the same actual DBI work and output fields.
local({
  state = new.env(parent = emptyenv()); state[["queries"]] = 0L
  state[["observe"]] = function() state[["queries"]] = state[["queries"]]+1L
  key = ".dbreg_null_query_count__"
  stopifnot(!exists(key,.GlobalEnv,inherits=FALSE))
  assign(key,state,.GlobalEnv)
  invisible(capture.output(trace("dbGetQuery",where=asNamespace("DBI"),
    tracer=quote(get(".dbreg_null_query_count__",envir=.GlobalEnv)[["observe"]]()),print=FALSE)))
  on.exit({invisible(capture.output(untrace("dbGetQuery",where=asNamespace("DBI"))));rm(list=key,envir=.GlobalEnv)})
  d = data.frame(x=seq_len(60)/60,y=sin(seq_len(60)))
  state[["queries"]]=0L
  implicit = dbglm(y~x,data=d)
  implicit_queries = state[["queries"]]
  state[["queries"]]=0L
  explicit = dbglm(y~x,data=d,null_model=FALSE)
  expect_identical(state[["queries"]],implicit_queries)
  expect_null(implicit[["null_diagnostics"]])
  expect_null(explicit[["null_diagnostics"]])
  expect_identical(names(gof(implicit)),c("logLik","deviance","AIC","BIC"))
  expect_identical(names(generics::glance(implicit)),names(generics::glance(explicit)))
  expect_identical(capture.output(print(implicit)),capture.output(print(explicit)))
  state[["queries"]] = 0L
  intercept_default = dbglm(y~1,data=d)
  intercept_queries = state[["queries"]]
  state[["queries"]] = 0L
  requested = dbglm(y~1,data=d,null_model=TRUE)
  expect_identical(state[["queries"]], intercept_queries)
  expect_identical(requested[["null_diagnostics"]][["fit"]][["evaluations"]],0L)
  after_fit = state[["queries"]]
  invisible(gof(requested)); invisible(generics::glance(requested)); invisible(capture.output(print(requested)))
  close(requested)
  expect_identical(state[["queries"]],after_fit)
  close(implicit);close(explicit);close(intercept_default)
  state[["queries"]] = 0L
  expect_error(dbglm(y~x,data=d,null_model=NA), "null_model.*TRUE or FALSE")
  expect_identical(state[["queries"]], 0L)
})

# Null sample is the already-pruned full-model population, never a new formula.
local({
  set.seed(9211)
  d = expand.grid(g=1:6,h=1:4,replicate=1:8)
  d[["id"]]=seq_len(nrow(d));d[["x"]]=rnorm(nrow(d));d[["w"]]=runif(nrow(d),.5,1.5)
  d[["o"]]=.1*sin(d[["id"]]);d[["cl"]]=d[["h"]]
  d[["use"]]=TRUE;d[["y"]]=.2*d[["x"]]+sin(d[["g"]])+d[["o"]]+rnorm(nrow(d))
  d[["x"]][1]=NA;d[["cl"]][2]=NA;d[["w"]][3]=0;d[["use"]][4]=FALSE
  extra=d[5,,drop=FALSE];extra[["id"]]=nrow(d)+1L;extra[["g"]]=99L;d=rbind(d,extra)
  keep=complete.cases(d[c("x","cl","w","o","y")]) & d[["w"]]>0 & d[["use"]] & d[["g"]]!=99L
  fit=dbglm(y~x|g+h,data=d,family=stats::gaussian(),weights=~w,offset=~o,vcov=~cl,
    subset=~use,fixef_rm="singletons",null_model=TRUE,retain="inference",row_id="id")
  expected=stats::lm(y~1,data=d[keep,],weights=w,offset=o)
  n=sum(keep);rss=sum(d[["w"]][keep]*stats::residuals(expected)^2)
  ll=-n/2*(log(2*pi)+log(rss/n)+1)+sum(log(d[["w"]][keep]))/2-1
  expect_equal(fit[["null_diagnostics"]][["nobs"]],n)
  expect_equal(fit[["null_diagnostics"]][["reported"]][["logLik"]],ll,tolerance=1e-8)
  record=fit[["null_diagnostics"]]
  changed=summary(fit,vcov="hc1")
  expect_identical(changed[["null_diagnostics"]],record)
  close(fit)
  # A reported zero denominator and nonfinite full MLE are independent statuses.
  value=dbreg:::dbglm_null_likelihood(-2,-1,TRUE,NULL,TRUE)
  expect_equal(value[["ratio"]],-1,tolerance=0)
  expect_identical(dbreg:::dbglm_null_likelihood(-2,0,TRUE,NULL,TRUE)[["ratio_reason"]],"zero_null_likelihood")
  expect_identical(dbreg:::dbglm_null_likelihood(Inf,-1,TRUE,NULL,TRUE)[["ratio_reason"]],"nonfinite_full_likelihood")
})

# A failed optional computation is data about availability, not a failed full fit.
local({
  conn=DBI::dbConnect(duckdb::duckdb(config=list(threads="1",temp_directory="")))
  on.exit(DBI::dbDisconnect(conn,shutdown=TRUE))
  set.seed(9223)
  d=data.frame(x=rnorm(120),o=rep(c(-.1,.1),60),w=rep(c(.5,1.5,2.5),40))
  d[["y"]]=rep(0:1,60)
  full=dbglm(y~x,conn=conn,data=d,family=stats::binomial(),weights=~w,offset=~o)
  before_coef=coef(full);before_ll=as.numeric(logLik(full))
  source=data.frame("__dbglm_y"=d[["y"]],"__dbglm_a"=d[["w"]]/full[["weight_scale"]],
    "__dbglm_o"=d[["o"]],check.names=FALSE)
  DBI::dbWriteTable(conn,"null_budget_fixture",source)
  null=dbreg:::dbglm_null_diagnostics
  environment(null)=new.env(parent=environment(null))
  environment(null)[["dbglm_irls"]]=function(...) stop("GLM did not converge within the iteration budget.")
  expect_warning({record = null(full,conn,"null_budget_fixture",full[["family"]],full[["control"]],TRUE,
    full[["weight_scale"]],0,fractional_weights=TRUE)},"Requested null model unavailable")
  expect_false(record[["fit"]][["available"]])
  expect_identical(record[["fit"]][["reason"]],"budget_exceeded")
  expect_true(is.na(record[["fit"]][["evaluations"]]))
  expect_identical(coef(full),before_coef)
  expect_identical(as.numeric(logLik(full)),before_ll)
  expect_null(full[["null_diagnostics"]])
  expect_true(DBI::dbIsValid(conn))
  DBI::dbRemoveTable(conn,"null_budget_fixture")
  close(full)
  separated=data.frame(x=rep(c(-1,1),each=12),y=rep(0:1,each=12))
  expect_warning({boundary = dbglm(y~x,conn=conn,data=separated,family=stats::binomial(),
    separation="prediction",null_model=TRUE)},"boundary probabilities")
  expect_identical(boundary[["null_diagnostics"]][["fit"]][["reason"]],"nonfinite_full_mle")
  expect_false(boundary[["null_diagnostics"]][["weighted"]][["ratio_available"]])
  expect_equal(unname(gof(boundary)[["null.available"]]), 0)
  expect_true(is.na(gof(boundary)[["pseudo.r.squared"]]))
  expect_identical(generics::glance(boundary)[["null.reason"]], "nonfinite_full_mle")
  expect_true(any(grepl("nonfinite_full_mle", capture.output(print(boundary)), fixed = TRUE)))
  close(boundary)
})

# The global null drops FE, and optional retention never changes its population.
local({
  set.seed(9231)
  d = expand.grid(g = 1:5, replicate = 1:24)
  d[["id"]] = seq_len(nrow(d))
  d[["x"]] = sin(d[["id"]])
  d[["w"]] = rep(c(.5, 1.5, 2.5), length.out = nrow(d))
  d[["o"]] = .1*cos(d[["id"]])
  d[["y"]] = d[["g"]]/2+d[["o"]]+rnorm(nrow(d))
  fit = dbglm(y ~ 1 | g, data = d, weights = ~w, offset = ~o,
    fixef_rm = "none", null_model = TRUE)
  reference = stats::lm(y ~ 1, data = d, weights = w, offset = o)
  rss = sum(d[["w"]]*stats::residuals(reference)^2)
  expected = -nrow(d)/2*(log(2*pi)+log(rss/nrow(d))+1)+sum(log(d[["w"]]))/2-1
  expect_equal(fit[["null_diagnostics"]][["reported"]][["logLik"]], expected, tolerance = 1e-8)
  expect_true(abs(fit[["null_diagnostics"]][["reported"]][["logLik"]]-as.numeric(logLik(fit))) > 1)
  expect_true(fit[["null_diagnostics"]][["constant"]])
  close(fit)
  baseline = dbglm(y ~ x | g, data = d, weights = ~w, offset = ~o,
    fixef_rm = "none", null_model = TRUE)
  record = baseline[["null_diagnostics"]]
  for (policy in c("prediction", "inference", "all")) {
    retained = dbglm(y ~ x | g, data = d, weights = ~w, offset = ~o,
      fixef_rm = "none", null_model = TRUE, retain = policy, row_id = "id")
    expect_equal(retained[["null_diagnostics"]], record, tolerance = 1e-8)
    restored = unserialize(serialize(retained, NULL))
    expect_identical(restored[["null_diagnostics"]], retained[["null_diagnostics"]])
    close(retained)
    expect_identical(restored[["null_diagnostics"]], retained[["null_diagnostics"]])
    close(restored)
  }
  close(baseline)
  # Each supported binary link also has an analytic no-offset constant.
  d[["y"]] = as.numeric((d[["id"]] %% 5) < 2)
  mu = sum(d[["w"]]*d[["y"]])/sum(d[["w"]])
  expected = sum(round(d[["w"]])*stats::dbinom(d[["y"]], 1, mu, log = TRUE))
  weighted = sum(d[["w"]]*stats::dbinom(d[["y"]], 1, mu, log = TRUE))
  for (link in c("logit", "probit", "cloglog")) {
    binary = dbglm(y ~ x, data = d, family = stats::binomial(link), weights = ~w, null_model = TRUE)
    expect_identical(binary[["null_diagnostics"]][["fit"]][["method"]], "binomial_constant")
    expect_equal(binary[["null_diagnostics"]][["reported"]][["logLik"]], expected, tolerance = 1e-8)
    expect_equal(binary[["null_diagnostics"]][["weighted"]][["logLik"]], weighted, tolerance = 1e-8)
    close(binary)
  }
  # The shared Gaussian primitive retains scale in logarithms.
  value = dbreg:::dbglm_gaussian_loglik(10, 1e200, 1e200, -5)
  expect_true(is.finite(value))
  expect_equal(value, -5*(log(2*pi)+log(1e200/10)+log(1e200)+1)-2.5, tolerance = 0)
})
