local({
  original_solver = getOption("dbreg.ap_solver")
  on.exit(options(dbreg.ap_solver = original_solver), add = TRUE)
  set.seed(71)
  d = data.frame(id = 1:360, x = rnorm(360), z = rnorm(360), g = rep(1:12, each = 30),
    h = sample(1:6, 360, TRUE), cl = rep(1:18, each = 20), w = rep(c(.2, .7, 1), 120), o = 0)
  d[["y"]] = .5*d[["x"]] - .2*d[["z"]] + sin(d[["g"]]) + rnorm(360)
  d[["b"]] = stats::rbinom(360, 1L, stats::plogis(.2*d[["x"]]))
  d[["y"]][2L] = NA_real_
  d[["b"]][2L] = NA_real_
  d[["w"]][3L] = 0
  d[["cl"]][4L] = NA_integer_
  models = list()
  on.exit(for (fit in models) close(fit), add = TRUE)
  fit = dbreg(y ~ x + z | g + h, data = d, weights = "w", strategy = "demean", retain = "all")
  models[[1L]] = fit
  reference = fixest::feols(y ~ x + z | g + h, d, weights = ~w, fixef.rm = "none", notes = FALSE)
  expect_identical(obs(fit), fixest::obs(reference))
  expect_identical(case.names(fit), fixest::obs(reference))
  expect_equal(predict(fit), as.numeric(predict(reference)), tolerance = 1e-7)
  expect_equal(predict(fit, sample = "original"), as.numeric(predict(reference, sample = "original")), tolerance = 1e-7)
  expect_equal(residuals(fit, na.rm = FALSE), as.numeric(residuals(reference, na.rm = FALSE)), tolerance = 1e-7)
  expect_equal(df.residual(fit), df.residual(reference))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)), tolerance = 1e-8)
  before = list(coef = coef(fit), vcov = vcov(fit), nobs = nobs(fit))
  for (G in c("min", "conventional")) for (K in c("full", "none", "nonnested")) {
    policy = fixest::ssc(K.fixef = K, G.df = G)
    actual = vcov(fit, vcov = ~g+h, ssc = policy)
    expected = vcov(reference, vcov = ~g+h, ssc = policy, attr = TRUE)
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-7)
    expect_equal(attr(actual, "df.K"), attr(expected, "df.K"))
    expect_equal(attr(actual, "df.t"), attr(expected, "df.t"))
  }
  for (policy in list(fixest::ssc(K.adj = FALSE, G.adj = FALSE), fixest::ssc(t.df = "conventional"), fixest::ssc(t.df = 7))) {
    actual = summary(fit, vcov = ~g+h, ssc = policy)
    expected = summary(reference, vcov = ~g+h, ssc = policy)
    expect_equal(as.numeric(actual[["coeftable"]]), as.numeric(expected[["coeftable"]]), tolerance = 1e-6)
    expect_equal(as.numeric(confint(actual)), as.numeric(as.matrix(confint(expected))), tolerance = 1e-6)
  }
  expect_identical(list(coef = coef(fit), vcov = vcov(fit), nobs = nobs(fit)), before)
  invisible(capture.output(print(summary(fit, vcov = ~g+h, ssc = fixest::ssc()))))
  conn = dbreg:::fitted_resource_connection(fit[["resource_state"]])
  expect_equal(length(DBI::dbListTables(conn)), 1L)
  data = d[1:8, c("id", "x", "z", "g", "h")]
  point = predict(fit, data)
  copy = fit
  close(fit)
  expect_false(DBI::dbIsValid(conn))
  expect_error(obs(copy), "released")
  expect_error(predict(copy), "released")
  expect_error(vcov(copy, vcov = "hc1"), "released")
  expect_equal(predict(copy, data), point)
  expect_equal(as.numeric(logLik(copy)), as.numeric(logLik(reference)), tolerance = 1e-8)

  # Retain each projection path's original keys, including forced iterative AP.
  for (method in c("moments", "iterative")) {
    options(dbreg.ap_solver = method)
    fit = dbreg(y ~ x + z | g + h, data = d, weights = "w", strategy = "demean", retain = "all")
    models[[length(models)+1L]] = fit
    expect_identical(obs(fit), fixest::obs(reference))
    expect_equal(predict(fit), as.numeric(predict(reference)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit, vcov = ~g+h, ssc = fixest::ssc())), as.numeric(vcov(reference, vcov = ~g+h)), tolerance = 1e-6)
    close(fit)
  }
  options(dbreg.ap_solver = NULL)

  # A fit-time cluster omission becomes part of the immutable fit sample.
  fit = dbreg(y ~ x + z | g + h, data = d, weights = "w", vcov = ~cl,
    strategy = "demean", retain = "all")
  models[[length(models)+1L]] = fit
  ref_cluster = fixest::feols(y ~ x + z | g + h, d, weights = ~w, vcov = ~cl, fixef.rm = "none", notes = FALSE)
  expect_identical(obs(fit), fixest::obs(ref_cluster))
  expect_equal(as.numeric(vcov(fit, vcov = ~cl, ssc = fixest::ssc())), as.numeric(vcov(ref_cluster)), tolerance = 1e-6)
  expect_equal(as.numeric(vcov(fit, ssc = fixest::ssc())), as.numeric(vcov(ref_cluster)), tolerance = 1e-6)
  close(fit)

  fit = dbglm(b ~ x + z | g + h, data = d, weights = "w", offset = "o",
    family = stats::binomial("logit"), fixef_rm = "none", retain = "all")
  models[[length(models)+1L]] = fit
  ref_glm = fixest::feglm(b ~ x + z | g + h, d, weights = ~w, offset = ~o,
    family = stats::binomial("logit"), fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
  expect_identical(obs(fit), fixest::obs(ref_glm))
  expect_equal(predict(fit, type = "response"), as.numeric(predict(ref_glm)), tolerance = 1e-6)
  expect_equal(residuals(fit, type = "working"), as.numeric(residuals(ref_glm, type = "working")), tolerance = 1e-6)
  expect_equal(residuals(fit, na.rm = FALSE), as.numeric(residuals(ref_glm, na.rm = FALSE)), tolerance = 1e-6)
  expect_error(vcov(fit, vcov = ~g+h), "one-way GLM")
  close(fit)

  # One explicitly nested FE does not imply every other FE is nonnested.
  nested = d
  nested[["h"]] = 2L*nested[["g"]] + as.integer((seq_len(nrow(nested))-1L) %% 30L >= 15L)
  nested[["copyg"]] = nested[["g"]]
  fit = dbglm(b ~ x + z | g + h, data = nested, weights = "w", vcov = ~g,
    family = stats::binomial("logit"), fixef_rm = "none", retain = "all")
  models[[length(models)+1L]] = fit
  reference = fixest::feglm(b ~ x + z | g + h, nested, weights = ~w, vcov = ~g,
    family = stats::binomial("logit"), fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
  expect_identical(fit[["fe_nested"]], c(TRUE, FALSE))
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-6)
  expect_equal(fit[["n_params_vcov"]], attr(vcov(reference, attr = TRUE), "df.K"))
  expect_equal(as.numeric(vcov(fit, vcov = ~copyg)), as.numeric(vcov(reference, vcov = ~copyg)), tolerance = 1e-6)
  expect_equal(attr(vcov(fit, vcov = ~copyg), "df.K"), attr(vcov(reference, vcov = ~copyg, attr = TRUE), "df.K"))
  close(fit)

  V = matrix(c(1, 2, 2, 1), 2L)
  expect_warning({ repaired = dbreg:::inference_repair_vcov(V) }, "not positive definite")
  expect_equal(as.numeric(repaired), rep(1.5, 4L), tolerance = 1e-14)
  expect_warning({ unchanged = dbreg:::inference_repair_vcov(V, FALSE) }, "not positive definite")
  expect_identical(unchanged, V)
  expect_error(dbreg:::inference_repair_vcov(matrix(Inf, 1L)), "nonfinite")
})
