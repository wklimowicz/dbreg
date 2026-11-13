local({
  set.seed(11154)
  n = 120L
  d = data.frame(id = sprintf("r%03d", seq_len(n)), g = rep(1:4, each = n / 4),
    x = rnorm(n), z = rnorm(n), w = runif(n, .25, 4), o = rnorm(n))
  d[["y"]] = .4 * d[["x"]] + cos(d[["g"]]) + d[["o"]] + rnorm(n)
  d[["x"]][5L] = NA
  d[["w"]][7L] = 0
  keep = !is.na(d[["x"]]) & d[["w"]] > 0
  ref = lm(y ~ x + z, data = d[keep, ], weights = w, offset = o)
  fit = dbreg(y ~ x + z, data = d, weights = ~w, offset = ~o,
    retain = "all", strategy = "moments")
  expect_equal(fitted(fit), unname(fitted(ref)), tolerance = 1e-10)
  expect_identical(case.names(fit), which(keep))
  expect_equal(residuals(fit), unname(residuals(ref)), tolerance = 1e-10)
  expect_equal(residuals(fit, "pearson"),
    unname(residuals(ref)) * sqrt(d[["w"]][keep]), tolerance = 1e-10)
  expect_identical(which(is.na(residuals(fit, na.rm = FALSE))), which(!keep))
  close(fit)
  expect_error(fitted(fit), "no live database session")

  d[["h"]] = rep(1:5, length.out = n)
  for (groups in list(y ~ x + z | g, y ~ x + z | g + h)) {
    two = length(all.vars(stats::formula(Formula::Formula(groups), lhs = 0, rhs = 2))) == 2L
    formula = if (two) y ~ x + z + factor(g) + factor(h) else y ~ x + z + factor(g)
    ref = lm(formula, data = d[keep, ], weights = w, offset = o)
    fit = dbreg(groups, data = d, weights = ~w, offset = ~o,
      retain = "inference", strategy = "demean")
    expect_equal(fitted(fit), unname(fitted(ref)), tolerance = 1e-10)
    expect_equal(residuals(fit), unname(residuals(ref)), tolerance = 1e-10)
    close(fit)
  }

  varying = dbreg(y ~ x | g[z], data = d, weights = ~w, offset = ~o,
    retain = "all", strategy = "demean")
  ref = lm(y ~ x + factor(g)*z, data = d[keep, ], weights = w, offset = o)
  expect_equal(fitted(varying), unname(fitted(ref)), tolerance = 1e-10)
  close(varying)

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "retained_ols", d)
  fit = dbreg(y ~ x + z | g, conn = conn, table = "retained_ols", row_id = "id",
    weights = ~w, offset = ~o, retain = "all", strategy = "mundlak")
  selected = d[keep, ]
  selected[["mx"]] = ave(selected[["w"]] * selected[["x"]], selected[["g"]], FUN = sum) /
    ave(selected[["w"]], selected[["g"]], FUN = sum)
  selected[["mz"]] = ave(selected[["w"]] * selected[["z"]], selected[["g"]], FUN = sum) /
    ave(selected[["w"]], selected[["g"]], FUN = sum)
  ref = lm(y ~ x + z + mx + mz, data = selected, weights = w, offset = o)
  values = fitted(fit)
  close(fit)
  result = as.data.frame(values)
  expect_equal(result[[".fitted"]], unname(fitted(ref)), tolerance = 1e-10)
  expect_identical(result[["id"]], d[["id"]][keep])
  close(values)
  expect_identical(DBI::dbListTables(conn), "retained_ols")
})
