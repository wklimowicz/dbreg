library(fixest)

local({
  set.seed(183)
  d = data.frame(x = rnorm(390), period = rep(-6:6, 30), g = rep(1:10, each = 39))
  d[["y"]] = 0.3 * d[["x"]] + 0.1 * d[["period"]] + rnorm(390)
  d[["b"]] = stats::rbinom(390, 1, stats::plogis(0.2 * d[["x"]] + 0.1 * d[["period"]]))
  base_yr = -2
  i = fixest::i
  fml = y ~ x + i(period, ref = base_yr, keep = -5:5) | g
  reference = fixest::feols(fml, data = d, fixef.rm = "none", notes = FALSE)
  linear = dbreg::dbreg(fml, data = d, strategy = "demean")
  gaussian = dbreg::dbglm(fml, data = d, family = stats::gaussian())
  expect_equal(stats::coef(linear), stats::coef(reference), tolerance = 1e-9)
  expect_equal(stats::coef(gaussian), stats::coef(reference), tolerance = 1e-9)
  fml = b ~ x + i(period, ref = base_yr, keep = -5:5) | g
  reference = fixest::feglm(fml, data = d, family = stats::binomial(), fixef.rm = "none", notes = FALSE)
  fit = dbreg::dbglm(fml, data = d, family = stats::binomial())
  expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-6)

  d[["f"]] = rep(c("alpha", "alpine", "Beta"), 130)
  cases = list(
    list(fml = y ~ x + i(f, ref = "@^al"), reference = y ~ x + i(f, ref = c("alpha", "alpine"))),
    list(fml = y ~ x + i(f, ref = "@!^al"), reference = y ~ x + i(f, ref = "Beta")),
    list(fml = y ~ x + i(period, -2), reference = y ~ x + i(period, ref = -2))
  )
  for (case in cases) {
    reference = fixest::feols(case[["reference"]], data = d, notes = FALSE)
    fit = dbreg::dbreg(case[["fml"]], data = d, strategy = "moments")
    gaussian = dbreg::dbglm(case[["fml"]], data = d, family = stats::gaussian())
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-9)
    expect_equal(stats::coef(gaussian), stats::coef(reference), tolerance = 1e-9)
  }
  # Compact formulas retain evaluated selectors, and R/SQL predictions share
  # the fitted i() coding even after the fitting environment has gone away.
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  d[["id"]] = seq_len(nrow(d))
  duckdb::duckdb_register(conn, "input", d)
  fml = y ~ x + i(period, ref = base_yr, keep = -5:5) | g
  fit = dbreg::dbreg(fml, conn = conn, table = "input", retain = "prediction", strategy = "demean")
  reference = fixest::feols(fml, data = d, fixef.rm = "none", notes = FALSE)
  expect_equal(stats::predict(fit, d), stats::predict(reference, d), tolerance = 1e-9)
  output = stats::predict(fit, conn = conn, table = "input", row_id = "id")
  on.exit(close(output), add = TRUE)
  values = as.data.frame(output)
  expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
    stats::predict(reference, d), tolerance = 1e-9)
  close(output)

  selector_state = new.env(parent = emptyenv())
  selector_state[["calls"]] = 0L
  selector = function() {
    selector_state[["calls"]] = selector_state[["calls"]] + 1L
    -2
  }
  fit = dbreg::dbreg(y ~ x + i(period, ref = selector()) | g,
    data = d, strategy = "demean", retain = "prediction")
  expect_equal(selector_state[["calls"]], 1L)
  reference = fixest::feols(y ~ x + i(period, ref = -2) | g,
    data = d, fixef.rm = "none", notes = FALSE)
  expect_equal(stats::predict(fit, d), stats::predict(reference, d), tolerance = 1e-9)
  expect_equal(selector_state[["calls"]], 1L)
  selector_state[["calls"]] = 0L
  d[["instrument"]] = d[["x"]] + stats::rnorm(nrow(d), sd = 0.1)
  fit = dbreg::dbivreg(y ~ i(period, ref = selector()) | x ~ instrument,
    data = d, strategy = "moments")
  expect_equal(selector_state[["calls"]], 1L)
  reference = fixest::feols(y ~ i(period, ref = -2) | x ~ instrument, data = d, notes = FALSE)
  names_reference = sub("^fit_", "", names(stats::coef(reference)))
  expect_equal(unname(stats::coef(fit)[names_reference]), unname(stats::coef(reference)), tolerance = 1e-8)

})
