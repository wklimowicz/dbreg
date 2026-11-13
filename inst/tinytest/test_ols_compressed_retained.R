local({
  set.seed(11155)
  d = expand.grid(g = 1:4, h = 1:3, x = c(-1, 1), z = c(-1, 1), rep = 1:6)
  n = nrow(d)
  d[["id"]] = sprintf("row-%03d", seq_len(n))
  d[["w"]] = runif(n, .2, 3)
  d[["o"]] = runif(n, -.4, .4)
  d[["y"]] = .5 * d[["x"]] + cos(d[["g"]]) * d[["z"]] + d[["o"]] + rnorm(n)
  d[["x"]][2L] = NA_real_
  d[["w"]][3L] = 0
  keep = !is.na(d[["x"]]) & d[["w"]] > 0
  cases = list(
    list(y ~ x + z + x:z | g, y ~ x + z + x:z + factor(g)),
    list(y ~ x | (g^h)[z], y ~ x + interaction(g, h)*z)
  )
  for (case in cases) {
    fit = dbreg(case[[1L]], data = d, weights = ~w, offset = ~o,
      strategy = "compress", retain = "all")
    ref = lm(case[[2L]], data = d[keep, ], weights = w, offset = o)
    expect_identical(fit[["strategy"]], "compress")
    expect_true(fit[["nobs"]] < sum(keep))
    expect_equal(fitted(fit), unname(fitted(ref)), tolerance = 1e-10)
    expect_equal(residuals(fit), unname(residuals(ref)), tolerance = 1e-10)
    expect_identical(case.names(fit), which(keep))
    expect_identical(which(is.na(residuals(fit, na.rm = FALSE))), which(!keep))
    close(fit)
  }
  constant = dbreg(y ~ 1, data = d, weights = ~w, retain = "all",
    strategy = "compress")
  expect_equal(fitted(constant), rep(weighted.mean(d[["y"]], d[["w"]]), sum(d[["w"]] > 0)))
  close(constant)
  auto = dbreg(y ~ x + z, data = d, weights = ~w, retain = "all", compress_ratio = 1)
  expect_identical(auto[["strategy"]], "compress")
  close(auto)

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "compressed_retain_source", d)
  fit = dbreg(y ~ x + z, conn = conn, table = "compressed_retain_source",
    row_id = "id", weights = ~w, strategy = "compress", retain = "all")
  expect_true(fit[["nobs"]] < sum(keep))
  result = fitted(fit)
  close(fit)
  expect_identical(as.data.frame(result)[["id"]], d[["id"]][keep])
  close(result)
  expect_identical(DBI::dbListTables(conn), "compressed_retain_source")
})
