local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11145)
  d = expand.grid(g = 1:5, h = 1:3, x = c(-1, 0, 1), rep = 1:12)
  n = nrow(d)
  d[["w"]] = sample(c(.5, 1.5, 2.5, 3.5), n, TRUE)
  d[["o"]] = d[["h"]] / 10
  d[["cl"]] = d[["g"]] + d[["h"]] %% 2L
  eta = .45 * d[["x"]] + cos(d[["g"]]) / 4 + sin(d[["h"]]) / 3 + d[["o"]]
  initial = DBI::dbListTables(conn)
  for (family in list(gaussian(), poisson(), binomial("logit"),
      binomial("probit"), binomial("cloglog"))) {
    mu = family[["linkinv"]](eta)
    d[["y"]] = switch(family[["family"]], gaussian = rnorm(n, mu),
      poisson = rpois(n, mu), binomial = rbinom(n, 1, mu))
    for (fml in list(y ~ x, y ~ x | g + h)) {
      kind = if (family[["link"]] %in% c("log", "logit")) ~cl else "hc1"
      fits = lapply(c("moments", "compress"), function(strategy) {
        dbglm(fml, data = d, conn = conn, family = family, weights = ~w,
          offset = ~o, vcov = kind, ssc = "full", null_model = TRUE,
          strategy = strategy)
      })
      row = fits[[1L]]
      cell = fits[[2L]]
      expect_equal(nobs(cell), n)
      expect_true(cell[["n_cells"]] < n / 4)
      expect_identical(cell[["strategy"]], "compress")
      expect_equal(coef(cell), coef(row), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(cell)), as.numeric(vcov(row)), tolerance = 1e-8)
      expect_equal(cell[["deviance"]], row[["deviance"]], tolerance = 1e-10)
      expect_equal(as.numeric(logLik(cell)), as.numeric(logLik(row)), tolerance = 1e-10)
      expect_equal(as.numeric(logLik(cell, type = "weighted")),
        as.numeric(logLik(row, type = "weighted")), tolerance = 1e-10)
      expect_equal(cell[["null_diagnostics"]][["reported"]][["logLik"]],
        row[["null_diagnostics"]][["reported"]][["logLik"]], tolerance = 1e-10)
      expect_identical(DBI::dbListTables(conn), initial)
    }
  }
  fit = dbglm(y ~ x, data = d, conn = conn, family = binomial(), compress_nmax = 1L)
  expect_identical(fit[["strategy"]], "moments")
  fit = dbglm(y ~ x, data = d, conn = conn, family = binomial())
  expect_identical(fit[["strategy"]], "compress")
  expect_error(dbglm(y ~ x, data = d, conn = conn, family = binomial(),
    strategy = "compress", compress_nmax = 1L), "exceeds the cell limit")
  expect_identical(DBI::dbListTables(conn), initial)
})
