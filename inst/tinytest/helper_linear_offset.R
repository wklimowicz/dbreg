dbreg_test_linear_offset = function(conn) {
  set.seed(541)
  n = 180L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rep(0:2, 60),
    g = rep(1:9, each = 20), h = rep(1:6, 30), o = rnorm(n), w = runif(n, .1, 1.7))
  d[["y"]] = 2 + .4*d[["x"]] - .3*d[["z"]] + d[["g"]]/5 + d[["o"]] + rnorm(n)
  d[["w"]][4] = 0
  d[["o"]][8] = NA_real_
  d[["x"]][12] = NA_real_
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  baseline = dbglm_validation_tables(conn)
  cases = list(list(fml = y ~ x + z, strategy = "moments"),
    list(fml = y ~ z | g, strategy = "compress"),
    list(fml = y ~ x + z | g, strategy = "demean"),
    list(fml = y ~ x + z | g + h, strategy = "demean"))
  for (case in cases) {
    fit = dbreg(case[["fml"]], conn = conn, table = table, weights = "w", offset = "o",
      retain = "all", row_id = "id", strategy = case[["strategy"]], diagnostics = TRUE)
    ref = fixest::feols(case[["fml"]], data = d, weights = ~w, offset = ~o, fixef.rm = "none",
      ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit, vcov = "hc1")), as.numeric(vcov(ref, vcov = "hetero")), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit, cluster = ~h)), as.numeric(vcov(ref, cluster = ~h, ssc = fixest::ssc(K.fixef = "full"))), tolerance = 1e-8)
    expect_equal(as.numeric(nobs(fit)), as.numeric(nobs(ref)))
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-8)
    ids = fixest::obs(ref)
    p = predict(fit, d)
    expect_equal(p, predict(ref, d), tolerance = 1e-7)
    output = predict(fit)
    values = as.data.frame(output)
    order = match(ids, values[["id"]])
    expect_equal(values[[output[["prediction_column"]]]][order], predict(ref), tolerance = 1e-7)
    close(output)
    r = residuals(fit)
    values = as.data.frame(r)
    expect_equal(values[[r[["value_column"]]]][match(ids, values[["id"]])], residuals(ref), tolerance = 1e-7)
    close(r)
    output = predict(fit, newdata = d, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])], p, tolerance = 1e-7)
    close(output)
    close(fit)
    expect_identical(dbglm_validation_tables(conn), baseline)
  }
}
