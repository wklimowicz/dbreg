# Reject statistical requests that cannot yet be honored, before touching data.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  expect_error(dbreg(y ~ x, conn = conn, table = "absent", offset = ~o + x), "offset.*must name")
  expect_error(dbreg(y ~ x, conn = conn, table = "absent", weigths = "w"), "Unused arguments:.*weigths")
  expect_error(dbivreg(y ~ x | z, conn = conn, table = "absent", weigths = "w"), "Unused arguments:.*weigths")
  expect_identical(DBI::dbListTables(conn), character())

  set.seed(531)
  d = data.frame(x = rnorm(120), g = rep(1:6, each = 20), cl = rep(1:12, each = 10))
  d[["y"]] = .5 * d[["x"]] + d[["g"]] + rnorm(120)
  m = dbreg(y ~ x | g, conn = conn, data = d, strategy = "compress")
  before = list(coef = coef(m), vcov = vcov(m), nobs = nobs(m))
  expect_error(vcov(m, cluster = ~cl), "retain")
  expect_error(vcov(m, ssc = list(K.adj = FALSE)), "retain")
  expect_error(vcov(m, FALSE, TRUE), "Unused arguments:.*unnamed")
  expect_error(coef(m, keep = "x"), "Unused arguments:.*keep")
  expect_error(summary(m, unknown = TRUE), "Unused arguments:.*unknown")
  expect_error(nobs(m, unknown = TRUE), "Unused arguments:.*unknown")
  expect_error(df.residual(m, unknown = TRUE), "Unused arguments:.*unknown")
  expect_error(logLik(m, unknown = TRUE), "Unused arguments:.*unknown")
  expect_error(fitted(m, se.fit = TRUE), "Unused arguments:.*se.fit")
  expect_error(residuals(m, se.fit = TRUE), "Unused arguments:.*se.fit")
  expect_error(case.names(m, unknown = TRUE), "Unused arguments:.*unknown")
  expect_error(confint(m, cluster = ~cl), "retain")
  expect_error(predict(m, d, type = "link"), "Unused arguments:.*type")
  expect_identical(list(coef = coef(m), vcov = vcov(m), nobs = nobs(m)), before)

  # Unsupported confidence levels reject before reading prediction inputs.
  plain = dbreg(y ~ x, data = d, conn = conn, strategy = "moments")
  for (value in list(0, 1, 2, NA_real_, Inf, c(.9, .95), "95%")) {
    for (interval in c("confidence", "prediction")) {
      expect_error(predict(plain, data.frame(), interval = interval, level = value), "level")
    }
  }
  dense = lm(y ~ x, data = d)
  for (interval in c("confidence", "prediction")) {
    expected = predict(dense, d, interval = interval, level = .9)
    actual = predict(plain, d, interval = interval, level = .9)
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = 1e-7)
  }
  expect_identical(DBI::dbListTables(conn), character())

  # Keep legacy aliases working and preserve their values and warnings.
  expect_warning({ with_fe = coef(m, fes = TRUE) }, "superseded")
  expect_equal(with_fe, coef(m, fe = TRUE))
  expect_warning({ ci = confint(m, fes = TRUE) }, "superseded")
  expect_equal(ci, confint(m, fe = TRUE))
  expect_warning({ full = vcov(m, fes = TRUE) }, "superseded")
  expect_equal(full, vcov(m, fe = TRUE))
  invisible(capture.output(expect_warning({ sql = dbreg(y ~ x, data = d, query_only = TRUE) }, "superseded")))
  expect_true(is.character(sql))

  d[["y"]] = rep(c(0, 1, 1, 0), 30)
  g = dbglm(y ~ x, conn = conn, data = d, family = binomial("logit"))
  expect_error(vcov(g, cluster = ~cl), "retain")
  expect_error(confint(g, cluster = ~cl), "retain")
  expect_error(confint(g, fe = TRUE), "Absorbed-effect confidence intervals")
  for (value in list(0, 1, NA_real_, Inf, c(.9, .95), "95%")) {
    expect_error(confint(g, level = value), "level")
  }
  expect_error(confint(g, parm = NA), "Unknown coefficient")
  expect_equal(tidy(g, conf.int = TRUE)[["conf.low"]], unname(confint(g)[, 1]))
  expect_error(predict(g, d, se.fit = TRUE), "Unused arguments:.*se.fit")
  expect_error(logLik(g, unknown = TRUE), "Unused arguments:.*unknown")
  expect_true(is.finite(as.numeric(logLik(g))))
  expect_identical(DBI::dbListTables(conn), character())
})
