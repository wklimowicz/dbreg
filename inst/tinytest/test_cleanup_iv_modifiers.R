local({
  set.seed(972)
  n = 360L
  d = data.frame(x = sample(0:2, n, TRUE), z = sample(0:2, n, TRUE),
    g = rep(1:12, each = 30), cl = sample(1:15, n, TRUE), o = rnorm(n),
    w = runif(n, 0.2, 1.2), keep = rep(c(TRUE, TRUE, FALSE), 120))
  d[["d"]] = d[["z"]] + sample(0:2, n, TRUE)
  d[["y"]] = 0.4 * d[["x"]] + 0.7 * d[["d"]] + sin(d[["g"]]) + d[["o"]] + rnorm(n)
  d[["o"]][2:3] = c(NA_real_, Inf)
  reference_data = d[d[["keep"]] & is.finite(d[["o"]]), ]
  adjusted_data = reference_data
  adjusted_data[["y"]] = adjusted_data[["y"]] - adjusted_data[["o"]]
  for (fixed in c(FALSE, TRUE)) {
    fml = if (fixed) y ~ x | g | d ~ z else y ~ x | d ~ z
    for (strategy in if (fixed) "demean" else c("moments", "compress")) for (kind in list("iid", "hc1", ~cl)) {
      fit = dbreg::dbivreg(fml, data = d, weights = ~w, offset = ~o, subset = ~keep,
        strategy = strategy, vcov = kind)
      reference = fixest::feols(fml, data = adjusted_data, weights = ~w,
        vcov = if (identical(kind, "hc1")) "hetero" else kind, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
      expected = stats::coef(reference)
      names(expected) = sub("^fit_", "", names(expected))
      names = names(stats::coef(fit))
      expect_equal(stats::coef(fit), expected[names], tolerance = 1e-7)
      V = stats::vcov(reference)
      dimnames(V) = lapply(dimnames(V), function(names) sub("^fit_", "", names))
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V[names, names]), tolerance = 1e-7)
      expect_equal(stats::nobs(fit), stats::nobs(reference))
      if (!fixed) {
        X = stats::model.matrix(~x + d, adjusted_data)
        Z = stats::model.matrix(~x + z, adjusted_data)
        w = adjusted_data[["w"]]
        P = crossprod(X, w * Z) %*% solve(crossprod(Z, w * Z))
        beta = solve(P %*% crossprod(Z, w * X), P %*% crossprod(Z, w * adjusted_data[["y"]]))
        expect_equal(unname(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
      }
      if (!fixed) expect_equal(stats::predict(fit, reference_data), as.numeric(stats::predict(reference, reference_data)) + reference_data[["o"]], tolerance = 1e-7)
    }
  }
})

local({
  set.seed(979)
  d = data.frame(x = rnorm(240), z = rnorm(240), g = rep(1:12, each = 20))
  d[["d"]] = d[["z"]] + rnorm(240)
  d[["offset value"]] = rnorm(240)
  d[["case weight"]] = runif(240, 0.3, 1.1)
  d[["keep row"]] = rep(c(TRUE, TRUE, FALSE), 80)
  d[["outcome value"]] = d[["x"]] + d[["d"]] + d[["offset value"]] + rnorm(240)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "iv_input", d)
  fit = dbreg::dbivreg(`outcome value` ~ x | g | d ~ z, conn = conn, table = "iv_input",
    weights = ~`case weight`, offset = ~`offset value`, subset = ~`keep row`)
  r = d[d[["keep row"]], ]
  r[["outcome value"]] = r[["outcome value"]] - r[["offset value"]]
  r[["ref_y"]] = r[["outcome value"]]
  r[["ref_weight"]] = r[["case weight"]]
  reference = fixest::feols(ref_y ~ x | g | d ~ z, data = r, weights = ~ref_weight, fixef.rm = "none", notes = FALSE)
  expected = stats::coef(reference)
  names(expected) = sub("^fit_", "", names(expected))
  expect_equal(stats::coef(fit), expected[names(stats::coef(fit))], tolerance = 1e-8)
  expect_identical(DBI::dbListTables(conn), "iv_input")
  expect_true(DBI::dbIsValid(conn))
  simple = dbreg::dbreg(`outcome value` ~ 1, data = d, offset = ~`offset value`, strategy = "moments")
  reference = stats::lm(`outcome value` ~ offset(`offset value`), data = d)
  expect_equal(stats::coef(simple), stats::coef(reference), tolerance = 1e-8)
  expect_equal(stats::predict(simple, d), as.numeric(stats::predict(reference, d)), tolerance = 1e-8)
})
