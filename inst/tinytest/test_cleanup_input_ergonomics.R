local({
  set.seed(945)
  d = data.frame(x = rnorm(240), z = rnorm(240), w = runif(240, 0.2, 1.2),
    o = rnorm(240, sd = 0.1), g = rep(1:12, each = 20))
  d[["d"]] = d[["z"]] + rnorm(240)
  d[["y"]] = 0.5 * d[["x"]] + 0.7 * d[["d"]] + d[["o"]] + rnorm(240)
  reference = stats::lm(y ~ x + d + offset(o), data = d, weights = w)
  for (family in list(NULL, "gaussian", stats::gaussian, stats::gaussian())) {
    args = list(fml = y ~ x + d, data = d, weights = ~w, offset = ~o)
    if (!is.null(family)) args[["family"]] = family
    fit = do.call(dbreg::dbglm, args)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_identical(fit[["family"]][["family"]], "gaussian")
    expect_identical(fit[["resource_state"]], NULL)
  }
  for (kind in c("binomial", "poisson")) {
    fun = getExportedValue("stats", kind)
    dd = d
    eta = 0.3 * dd[["x"]] + dd[["o"]]
    dd[["y"]] = if (kind == "binomial") stats::rbinom(nrow(dd), 1, stats::plogis(eta)) else
      stats::rpois(nrow(dd), exp(eta))
    glm_reference = stats::glm(y ~ x + offset(o), data = dd, family = fun())
    for (family in list(kind, fun, fun())) {
      fit = dbreg::dbglm(y ~ x, data = dd, family = family, offset = ~o)
      expect_equal(stats::coef(fit), stats::coef(glm_reference), tolerance = 1e-6)
      expect_identical(fit[["family"]][["family"]], kind)
    }
  }
  expect_error(dbreg::dbglm(y ~ x, data = d, family = "unknown"), "Supported family")
  expect_error(dbreg::dbreg(y ~ x, data = d, weights = ~w + x), "must name one column")
  expect_error(dbreg::dbglm(y ~ x, data = d, offset = ~o + x), "must name one column")
  if (requireNamespace("dbplyr", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
    DBI::dbWriteTable(conn, "ergonomic_input", d)
    lazy = dplyr::tbl(conn, "ergonomic_input")
    for (source in list(list(data = lazy), list(table = lazy))) {
      fit = do.call(dbreg::dbreg, c(list(fml = y ~ x + d, weights = ~w, offset = ~o), source))
      expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
      fit = do.call(dbreg::dbglm, c(list(fml = y ~ x + d, weights = ~w, offset = ~o), source))
      expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
      fit = do.call(dbreg::dbivreg, c(list(fml = y ~ x | d ~ z, weights = ~w), source))
      iv_reference = fixest::feols(y ~ x | d ~ z, data = d, weights = ~w, fixef.rm = "none", notes = FALSE)
      expected = stats::coef(iv_reference)
      names(expected) = sub("^fit_", "", names(expected))
      expect_equal(stats::coef(fit), expected[names(stats::coef(fit))], tolerance = 1e-8)
      bins = do.call(dbreg::dbbinsreg, c(list(fml = y ~ x, nbins = 8, ci = FALSE, plot = FALSE), source))
      expect_true(inherits(bins, "dbbinsreg"))
      expect_true(DBI::dbIsValid(conn))
    }
    ignored = dplyr::tbl(conn, dbplyr::sql("SELECT * FROM ergonomic_input WHERE x > 0"))
    fit = dbreg::dbreg(y ~ x, data = ignored, table = lazy)
    expect_equal(stats::coef(fit), stats::coef(stats::lm(y ~ x, d)), tolerance = 1e-8)
    fit = dbreg::dbreg(y ~ x, data = data.frame(unrelated = 1), table = lazy, weights = ~w)
    expect_equal(stats::coef(fit), stats::coef(stats::lm(y ~ x, d, weights = w)), tolerance = 1e-8)
  }
})
