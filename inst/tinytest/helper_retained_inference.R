dbreg_test_retained_inference = function(conn) {
  before = dbglm_validation_tables(conn)
  resources = dbreg:::dbglm_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::dbglm_release_resources(resources), add = TRUE)
  original = options(dbreg.prediction_max_entries = 1e6)
  on.exit(options(original), add = TRUE)
  set.seed(394)
  n = 480L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:12, each = 40),
    h = sample(1:6, n, TRUE), cl = rep(1:24, each = 20), w = rep(c(.2, .7, 1), length.out = n), o = sin(seq_len(n))*.1)
  d[["y"]] = 1 + .4*d[["x"]] - .2*d[["z"]] + sin(d[["g"]]) + rnorm(n, sd = .3)
  d[["b"]] = stats::rbinom(n, 1L, stats::plogis(.3*d[["x"]] + .1*sin(d[["g"]])))
  d[["count"]] = stats::rpois(n, exp(.2*d[["x"]] + .1*cos(d[["g"]])))
  d[["x"]][2L] = NA_real_
  d[["w"]][3L] = 0
  d[["cl"]][4L] = NA_integer_
  d[c("y", "b", "count")][5L, ] = NA_real_
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  fits = list()
  on.exit(for (fit in fits) close(fit), add = TRUE)
  cases = list(list(fml = y ~ x + z, strategy = "moments"),
    list(fml = y ~ x + z | g, strategy = "demean"),
    list(fml = y ~ x + z | g + h, strategy = "demean"),
    list(fml = y ~ x + z | g + h, strategy = "compress"))
  for (family in list(stats::binomial("logit"), stats::binomial("cloglog"), stats::binomial("probit"), stats::poisson())) {
    response = if (family[["family"]] == "poisson") "count" else "b"
    cases[[length(cases)+1L]] = list(fml = stats::as.formula(paste0(response, "~x+z|g+h")), family = family)
  }
  close_outputs = function(output) {
    on.exit(close(output))
    actual = as.data.frame(output)
    actual[order(actual[["id"]]), , drop = FALSE]
  }
  for (i in seq_along(cases)) {
    case = cases[[i]]
    glm = !is.null(case[["family"]])
    label = if (glm) case[["family"]][["link"]] else paste(case[["strategy"]], i)
    message("[Retained inference] ", label)
    args = list(fml = case[["fml"]], conn = conn, table = table, weights = "w", row_id = "id", retain = "all")
    if (glm) {
      args = c(args, list(family = case[["family"]], offset = "o", fixef_rm = "none"))
      reference = fixest::feglm(case[["fml"]], d, weights = ~w, offset = ~o,
        family = case[["family"]], fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
    } else {
      args[["strategy"]] = case[["strategy"]]
      reference = fixest::feols(case[["fml"]], d, weights = ~w, fixef.rm = "none", notes = FALSE)
    }
    fit = do.call(if (glm) dbglm else dbreg, args)
    fits[[i]] = fit
    state = fit[["fitted_state"]][["inference"]]
    expect_equal(state[["nobs"]], nobs(reference), info = label)
    expect_equal(state[["nobs_input"]], n, info = label)
    expect_true(as.numeric(object.size(fit)) < 512*1024, info = label)
    ids = close_outputs(dbreg::obs(fit))[["id"]]
    expect_equal(ids, fixest::obs(reference), info = label)
    prediction = if (glm) predict(fit, type = "link") else predict(fit)
    actual = close_outputs(prediction)
    expect_equal(actual[[prediction[["prediction_column"]]]], as.numeric(predict(reference, type = "link")), tolerance = 1e-6, info = label)
    original_prediction = if (glm) predict(fit, type = "response", sample = "original") else predict(fit, sample = "original")
    actual = close_outputs(original_prediction)
    expected = as.numeric(predict(reference, type = "response", sample = "original"))
    expect_identical(is.na(actual[[original_prediction[["prediction_column"]]]]), is.na(expected), info = label)
    expect_true(all(abs(actual[[original_prediction[["prediction_column"]]]] - expected) <= 1e-7 + 1e-6*abs(expected), na.rm = TRUE), info = label)
    for (kind in c("response", "pearson", "deviance", if (glm) "working")) {
      output = residuals(fit, type = kind)
      actual = close_outputs(output)
      expected = as.numeric(residuals(reference, type = kind))
      expect_true(all(abs(actual[[output[["value_column"]]]] - expected) <= 1e-7 + 1e-6*abs(expected)), info = paste(label, kind))
    }
    expect_true(abs(as.numeric(logLik(fit)) - as.numeric(logLik(reference))) <= 1e-7 + 1e-7*abs(as.numeric(logLik(reference))), info = label)
    expect_equal(attr(logLik(fit), "df"), attr(logLik(reference), "df"), info = label)
    requests = list("iid", "hc1", ~g, if (!glm) ~g+h)
    requests = requests[!vapply(requests, is.null, logical(1))]
    for (request in requests) {
      policy = fixest::ssc()
      V = vcov(fit, vcov = request, ssc = policy)
      expected = stats::vcov(reference, vcov = request, ssc = policy, attr = TRUE)
      expect_true(all(abs(V - expected) <= 1e-8 + 1e-5*abs(expected)), info = label)
      expect_equal(attr(V, "df.K"), attr(expected, "df.K"), info = label)
      expect_equal(attr(V, "df.t"), attr(expected, "df.t"), info = label)
    }
    selected = if (glm) ~g else ~g+h
    summary = summary(fit, vcov = selected, ssc = fixest::ssc())
    expected = summary(reference, vcov = selected)
    expect_true(all(abs(as.numeric(summary[["coeftable"]]) - as.numeric(expected[["coeftable"]])) <=
      1e-7 + 1e-5*abs(as.numeric(expected[["coeftable"]]))), info = label)
    expect_equal(as.numeric(confint(summary)), as.numeric(as.matrix(confint(expected))), tolerance = 1e-5, info = label)
    expect_equal(as.numeric(vcov(summary, ssc = list(K.adj = FALSE))),
      as.numeric(vcov(reference, vcov = selected, ssc = fixest::ssc(K.adj = FALSE))), tolerance = 1e-5, info = label)
    expect_error(vcov(fit, cluster = ~cl), "missing values on the fitted sample")
    expect_error(vcov(fit, cluster = ~absent), "absent from the retained source")
    expect_error(vcov(fit, ssc = list(K.exact = TRUE)), "Unsupported.*ssc")
    serialized = unserialize(serialize(fit, NULL))
    expect_error(predict(serialized), "live owning session")
    close(serialized)
  }
  message("[Retained inference] zero-column Poisson")
  fit = dbglm(count ~ 1 | g + h, conn = conn, table = table, weights = "w", offset = "o",
    family = stats::poisson(), fixef_rm = "none", row_id = "id", retain = "all")
  fits[[length(fits)+1L]] = fit
  reference = fixest::feglm(count ~ 1 | g + h, d, weights = ~w, offset = ~o,
    family = stats::poisson(), fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
  output = predict(fit)
  actual = close_outputs(output)
  expect_equal(actual[[output[["prediction_column"]]]], as.numeric(predict(reference, type = "link")), tolerance = 1e-6)
  expect_equal(dim(vcov(fit, vcov = "hc1")), c(0L, 0L))
  expect_equal(dim(vcov(fit, vcov = ~g)), c(0L, 0L))
  expect_equal(fit[["fitted_state"]][["inference"]][["K_model"]], reference[["nparams"]])
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)), tolerance = 1e-7)
  message("[Retained inference] all nested GLM fixed effects")
  nested = d
  nested[["h"]] = 2L*nested[["g"]] + as.integer((nested[["id"]]-1L) %% 40L >= 20L)
  nested[["copyg"]] = nested[["g"]]
  nested_table = dbreg:::prediction_upload_frame(conn, nested, resources)
  fit = dbglm(b ~ x + z | g + h, conn = conn, table = nested_table, weights = "w", vcov = ~g,
    family = stats::binomial("logit"), fixef_rm = "none", row_id = "id", retain = "all")
  fits[[length(fits)+1L]] = fit
  reference = fixest::feglm(b ~ x + z | g + h, nested, weights = ~w, vcov = ~g,
    family = stats::binomial("logit"), fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
  expect_identical(fit[["fe_nested"]], c(TRUE, FALSE))
  expect_true(all(abs(vcov(fit)-vcov(reference)) <= 1e-8 + 1e-5*abs(vcov(reference))))
  expect_equal(fit[["n_params_vcov"]], attr(vcov(reference, attr = TRUE), "df.K"))
  expect_true(all(abs(vcov(fit, vcov = ~copyg)-vcov(reference, vcov = ~copyg)) <= 1e-8 + 1e-5*abs(vcov(reference, vcov = ~copyg))))
  expect_equal(attr(vcov(fit, vcov = ~copyg), "df.K"), attr(vcov(reference, vcov = ~copyg, attr = TRUE), "df.K"))
  # No operation may silently reread a source replaced after fitting.
  dbglm_validation_drop(conn, gsub('["`]', "", table))
  for (i in c(1L, 3L, 8L)) {
    fit = fits[[i]]
    actual = close_outputs(predict(fit))
    expect_equal(nrow(actual), nobs(fit))
    expect_true(all(is.finite(vcov(fit, vcov = "hc1"))))
  }
  for (fit in fits) close(fit)
  for (fit in fits) expect_error(predict(fit), "released")
  dbreg:::dbglm_release_resources(resources)
  resources[["tables"]] = character()
  expect_identical(dbglm_validation_tables(conn), before)
  expect_true(DBI::dbIsValid(conn))
}
