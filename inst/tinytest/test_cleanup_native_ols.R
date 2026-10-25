local({
  set.seed(512)
  d = data.frame(x = rnorm(360), z = rnorm(360), f = rep(letters[1:3], 120), g = rep(1:12, each = 30))
  d[["y"]] = 0.3 * d[["x"]] + d[["g"]] / 10 + d[["z"]] * sin(d[["g"]]) + rnorm(360)
  cases = list(
    list(fml = y ~ 1, strategies = c("moments", "compress")),
    list(fml = y ~ x + factor(f), strategies = c("moments", "compress")),
    list(fml = y ~ 0 + factor(f), strategies = c("moments", "compress")),
    list(fml = y ~ 1 | g, strategies = c("demean", "compress")),
    list(fml = y ~ x + factor(f) | g, strategies = c("demean", "compress")),
    list(fml = y ~ x + i(f) | g, strategies = c("demean", "compress")),
    list(fml = y ~ x | g[z], strategies = c("demean", "compress")),
    list(fml = y ~ x + z | g[z], strategies = c("demean", "compress"))
  )
  expect_error(dbreg::dbreg(y ~ x | g[z], data = d, strategy = "demean", sql_only = TRUE),
    "SQL-only projection plans")
  expect_error(dbreg::dbreg(y ~ x | g[z], data = d, strategy = "mundlak"),
    "demean or compress")
  for (case in cases) {
    reference = fixest::feols(case[["fml"]], data = d, fixef.rm = "none", notes = FALSE)
    for (strategy in case[["strategies"]]) {
      fit = dbreg::dbreg(case[["fml"]], data = d, strategy = strategy)
      expect_true(inherits(fit, "dbreg") && !inherits(fit, "dbglm"))
      if (length(stats::coef(reference))) {
        expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
      } else expect_equal(length(stats::coef(fit)), 0L)
      expect_equal(stats::predict(fit, d), stats::predict(reference, d), tolerance = 1e-8)
      expect_equal(nrow(dbreg::glance(fit)), 1L)
      expect_true(all(is.finite(unlist(dbreg::glance(fit)[c("r.squared", "adj.r.squared", "rmse")]))))
    }
  }
})

local({
  set.seed(628)
  n = 480L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), z2 = rnorm(n),
    g = rep(1:16, each = 30), h = sample.int(7, n, TRUE), f = rep(letters[1:3], 160),
    w = rep(c(0.2, 0.7, 1.3), 160), cl = sample.int(20, n, TRUE))
  d[["y"]] = 0.4 * d[["x"]] + d[["g"]] / 8 + sin(d[["g"]]) * d[["z"]] +
    0.2 * d[["h"]] + rnorm(n)
  d[["x"]][2] = NA_real_
  d[["z"]][4] = NA_real_
  d[["y"]][6] = NA_real_
  d[["w"]][8] = 0
  formulas = list(y ~ x | g[z], y ~ x | g[z] + h,
    y ~ x + factor(f) | g[z, z2] + h, y ~ x | g[z, z2] + h[z])
  for (fml in formulas) for (strategy in c("demean", "compress")) {
    for (kind in list("iid", "hc1", ~g, ~cl)) {
      fit = dbreg::dbreg(fml, data = d, weights = "w", strategy = strategy, vcov = kind)
      reference = fixest::feols(fml, data = d, weights = ~w,
        vcov = if (identical(kind, "hc1")) "hetero" else kind,
        ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
      expect_equal(stats::nobs(fit), stats::nobs(reference))
      expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
      expect_equal(stats::df.residual(fit), stats::df.residual(reference))
    }
    fit = dbreg::dbreg(fml, data = d, weights = "w", strategy = strategy, retain = "all", vcov = ~cl)
    on.exit(close(fit), add = TRUE)
    reference = fixest::feols(fml, data = d, weights = ~w, fixef.rm = "none",
      ssc = fixest::ssc(K.fixef = "full"), fixef.tol = 1e-11, notes = FALSE)
    expect_equal(stats::predict(fit, d), stats::predict(reference, d), tolerance = 1e-7)
    expect_equal(stats::predict(fit), stats::predict(reference), tolerance = 1e-7)
    expect_equal(as.numeric(stats::vcov(fit, vcov = ~g)), as.numeric(stats::vcov(reference, vcov = ~g,
      ssc = fixest::ssc(K.fixef = "full"))), tolerance = 1e-7)
    conn = dbreg:::fitted_resource_connection(fit[["resource_state"]])
    output = stats::predict(fit, d, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
      stats::predict(reference, d), tolerance = 1e-7)
    close(output)
    close(fit)
  }
})

local({
  set.seed(724)
  d = data.frame(id = seq_len(240), y = rnorm(240), x = rnorm(240), z = rnorm(240))
  d[["category code"]] = rep(1:3, 80)
  d[["group id"]] = rep(1:12, each = 20)
  fml = y ~ x * factor(`category code`) | `group id`[z]
  reference = stats::lm(y ~ x * factor(`category code`) + factor(`group id`) +
    factor(`group id`):z, data = d)
  for (strategy in c("demean", "compress")) {
    fit = dbreg::dbreg(fml, data = d, strategy = strategy, retain = "prediction")
    expect_equal(stats::coef(fit), stats::coef(reference)[names(stats::coef(fit))], tolerance = 1e-8)
    expect_equal(stats::predict(fit, d[1:2, ]), as.numeric(stats::predict(reference, d[1:2, ])), tolerance = 1e-8)
    conn = DBI::dbConnect(duckdb::duckdb())
    output = stats::predict(fit, d, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
      as.numeric(stats::predict(reference, d)), tolerance = 1e-8)
    close(output)
    DBI::dbDisconnect(conn, shutdown = TRUE)
  }
})

local({
  set.seed(728)
  d = data.frame(id = 1:90, y = rnorm(90), day = as.Date("2020-01-01") + rep(1:3, 30))
  reference = stats::lm(y ~ factor(day), d)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbreg(y ~ factor(day), data = d, strategy = strategy)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-10)
    expect_equal(stats::predict(fit, d[1:2, ]), as.numeric(stats::predict(reference, d[1:2, ])), tolerance = 1e-10)
    output = stats::predict(fit, d, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
      as.numeric(stats::predict(reference, d)), tolerance = 1e-10)
    close(output)
  }
})

local({
  set.seed(741)
  d = data.frame(y = rnorm(480), x = rnorm(480), g = rep(1:16, each = 30),
    period = rep(1:8, 60), w = rep(c(0.3, 0.8, 1), 160))
  fml = y ~ x | g[period] + period
  reference = fixest::feols(fml, data = d, weights = ~w, fixef.rm = "none",
    fixef.tol = 1e-11, notes = FALSE)
  for (strategy in c("demean", "compress")) {
    fit = dbreg::dbreg(fml, data = d, weights = "w", strategy = strategy, retain = "all")
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-8)
    expect_equal(stats::predict(fit, d[1:2, ]), stats::predict(reference, d[1:2, ]), tolerance = 1e-8)
    expect_equal(stats::predict(fit), stats::predict(reference), tolerance = 1e-8)
    close(fit)
  }
})

local({
  set.seed(177)
  d = data.frame(id = 1:600, x = rnorm(600), g = rep(1:20, each = 30),
    z = 1e9 + rnorm(600), w = rep(c(0.2, 0.5, 1), 200))
  d[["y"]] = 0.4 * d[["x"]] + (d[["z"]] - 1e9) * sin(d[["g"]]) + rnorm(600)
  centered = d
  centered[["z"]] = centered[["z"]] - 1e9
  reference = fixest::feols(y ~ x | g[z], data = centered, weights = ~w,
    vcov = ~g, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
  for (strategy in c("demean", "compress")) {
    fit = dbreg::dbreg(y ~ x | g[z], data = d, weights = "w", strategy = strategy, vcov = ~g)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-8)
    expect_equal(stats::predict(fit, d), stats::predict(reference, centered), tolerance = 1e-7)
    if (strategy == "compress") {
      expect_equal(stats::predict(fit), stats::predict(fit, fit[["data"]]), tolerance = 1e-8)
    }
    retained = dbreg::dbreg(y ~ x | g[z], data = d, weights = "w", strategy = strategy,
      vcov = ~g, retain = "all")
    expect_equal(stats::predict(retained, d), stats::predict(reference, centered), tolerance = 1e-7)
    expect_equal(stats::predict(retained), stats::predict(reference, centered), tolerance = 1e-7)
    conn = dbreg:::fitted_resource_connection(retained[["resource_state"]])
    output = stats::predict(retained, d, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
      stats::predict(reference, centered), tolerance = 1e-7)
    close(output)
    close(retained)
  }
})
