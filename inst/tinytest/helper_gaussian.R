dbreg_test_gaussian = function(conn) {
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  models = outputs = list()
  on.exit({ for (output in outputs) close(output); for (model in models) close(model) }, add = TRUE)
  set.seed(567)
  n = 240L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:8, each = 30),
    h = rep(1:6, 40), trend = rnorm(n), o = rnorm(n), w = runif(n, .05, 1.5))
  d[["y"]] = 1+.4*d[["x"]]-.2*d[["z"]]+sin(d[["g"]])+d[["g"]]/8*d[["trend"]]+d[["o"]]+rnorm(n)
  d[["w"]][3] = 0
  d[["x"]][7] = NA_real_
  d[["o"]][11] = NA_real_
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  baseline = dbglm_validation_tables(conn)
  formulas = list(y ~ x + z, y ~ x + z | g, y ~ x + z | g+h, y ~ x+z | g[trend])
  for (fml in formulas) {
    message("[Gaussian] ", paste(deparse(fml), collapse = ""))
    model = dbglm(fml, conn = conn, table = table, family = stats::gaussian(), weights = "w", offset = "o",
      retain = "all", row_id = "id", fixef_rm = "none")
    models[[length(models)+1L]] = model
    ref = fixest::feglm(fml, d, family = "gaussian", weights = ~w, offset = ~o, fixef.rm = "none", notes = FALSE)
    expect_equal(coef(model), coef(ref), tolerance = 1e-8)
    expect_equal(as.numeric(logLik(model)), as.numeric(logLik(ref)), tolerance = 1e-8)
    expect_equal(attr(logLik(model), "df"), attr(logLik(ref), "df"))
    expect_equal(as.numeric(vcov(model)), as.numeric(vcov(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(model[["coeftable"]]), as.numeric(as.matrix(summary(ref)[["coeftable"]])), tolerance = 1e-7)
    expect_equal(as.numeric(confint(model)), as.numeric(as.matrix(confint(ref))), tolerance = 1e-7)
    expect_identical(model[["iterations"]], 1L)
    ids = fixest::obs(ref)
    for (request in list("iid", "hetero", ~g)) for (policy in list(fixest::ssc(), fixest::ssc(K.adj = FALSE), fixest::ssc(K.fixef = "none"))) {
      actual = summary(model, vcov = request, ssc = policy)
      expected = summary(ref, vcov = request, ssc = policy)
      expect_equal(as.numeric(vcov(actual)), as.numeric(vcov(expected)), tolerance = 1e-7)
      expect_equal(as.numeric(actual[["coeftable"]]), as.numeric(as.matrix(expected[["coeftable"]])), tolerance = 1e-6)
      expect_equal(as.numeric(confint(actual)), as.numeric(as.matrix(confint(expected))), tolerance = 1e-6)
    }
    for (type in c("response", "pearson", "deviance", "working")) {
      output = residuals(model, type = type)
      values = as.data.frame(output)
      expect_equal(values[[output[["value_column"]]]][match(ids, values[["id"]])],
        residuals(ref, type = type), tolerance = 1e-7)
      close(output)
    }
    point = predict(model, d, type = "response")
    expect_equal(point, predict(ref, d), tolerance = 1e-7)
    expect_equal(predict(model, d, type = "link"), point)
    output = predict(model, newdata = d, row_id = "id", type = "response")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])], point, tolerance = 1e-7)
    close(output)
    output = fitted(model)
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(ids, values[["id"]])], predict(ref), tolerance = 1e-7)
    close(output)
    expect_error(vcov(model, cluster = ~g+h), "one-way GLM")
    close(model)
    expect_equal(predict(model, d, type = "response"), point)
    expect_identical(dbglm_validation_tables(conn), baseline)
  }
  # The IID no-K policy uses RSS/(N-1), not RSS/N, in fixest OLS.
  ols = dbreg(y ~ x+z | g+h, conn = conn, table = table, weights = "w", offset = "o", retain = "all")
  models[[length(models)+1L]] = ols
  ref = fixest::feols(y ~ x+z | g+h, d, weights = ~w, offset = ~o, fixef.rm = "none", notes = FALSE)
  for (policy in list(fixest::ssc(K.adj = FALSE), fixest::ssc(K.adj = FALSE, K.fixef = "none"))) {
    expect_equal(as.numeric(vcov(ols, vcov = "iid", ssc = policy)), as.numeric(vcov(ref, vcov = "iid", ssc = policy)), tolerance = 1e-8)
  }
  close(ols)
  singleton = rbind(d, d[2, ])
  singleton[["id"]][n+1L] = n+1L
  singleton[["g"]][n+1L] = 99L
  pruning_table = dbreg:::prediction_upload_frame(conn, singleton, resources)
  fit = dbglm(y ~ x+z | g, conn = conn, table = pruning_table, family = stats::gaussian(),
    weights = "w", offset = "o", retain = "all", row_id = "id")
  models[[length(models)+1L]] = fit
  ref = fixest::feglm(y ~ x+z | g, singleton, family = "gaussian", weights = ~w, offset = ~o, notes = FALSE)
  expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
  expect_equal(nobs(fit), nobs(ref))
  expect_equal(fit[["nobs_pruned"]], 1)
  output = stats::case.names(fit)
  values = as.data.frame(output)
  expect_equal(sort(values[["id"]]), fixest::obs(ref))
  close(output)
  close(fit)
  # Uniform Gaussian weight scaling preserves all covariance types. A unit-
  # weight oracle avoids fixest's tiny-weight constant-outcome rejection.
  unit = d
  unit[["w"]] = 1
  reference = fixest::feglm(y ~ x+z | g, unit, family = "gaussian", offset = ~o, fixef.rm = "none", notes = FALSE)
  for (scale in c(1e-200, 1, 1e200)) {
    sql = paste("SELECT *,", dbreg:::dbglm_number(scale), "AS scale_weight FROM", table)
    name = basename(tempfile("dbreg_gaussian_scale_"))
    resources[["tables"]] = c(resources[["tables"]], name)
    DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, name, sql, c(names(d), "scale_weight")))
    fit = dbglm(y ~ x+z | g, conn = conn, table = name, family = stats::gaussian(), offset = "o", weights = "scale_weight", retain = "all", fixef_rm = "none")
    models[[length(models)+1L]] = fit
    expect_equal(coef(fit), coef(reference), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-7)
    for (request in list("iid", "hc1", ~g)) expect_equal(as.numeric(vcov(fit, vcov = request)),
      as.numeric(vcov(reference, vcov = if (identical(request, "hc1")) "hetero" else request)), tolerance = 1e-7)
    expect_true(is.finite(as.numeric(logLik(fit))))
    close(fit)
  }
}
