local({
  set.seed(1027)
  n = 720L
  d = data.frame(x = sample(0:3, n, TRUE), z = sample(0:2, n, TRUE),
    g = paste0("café '\"", sample(1:4, n, TRUE)),
    h = as.Date("2026-01-01") + sample(1:3, n, TRUE),
    c1 = sample(1:15, n, TRUE), c2 = sample(1:13, n, TRUE),
    w = stats::runif(n, 0.2, 1.4), o = stats::runif(n, -0.2, 0.2))
  d[["y"]] = 1 + 0.2 * d[["x"]] + 0.4 * d[["z"]] + d[["o"]] + stats::rnorm(n)
  for (formula in list(y ~ x + z, y ~ x + z | g + h, y ~ i(x, ref = 0) + z | g + h)) {
    reference = fixest::feols(formula, data = d, weights = ~w, offset = ~o,
      fixef.rm = "none", vcov = "hc1", ssc = fixest::ssc(K.fixef = "full"))
    for (mode in c("inference", "all")) {
      conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
      fit = dbreg(formula, conn = conn, data = d, strategy = "compress",
        weights = ~w, offset = ~o, retain = mode, vcov = "hc1")
      expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)), tolerance = 1e-7)
      expect_equal(stats::predict(fit), as.numeric(stats::predict(reference)), tolerance = 1e-7)
      expect_equal(stats::residuals(fit), as.numeric(stats::residuals(reference)), tolerance = 1e-7)
      for (kind in list("iid", "hc1", ~c1, ~c1 + c2)) {
        expect_equal(as.numeric(stats::vcov(fit, vcov = kind)),
          as.numeric(stats::vcov(reference, vcov = kind)), tolerance = 1e-7)
      }
      state = fit[["fitted_state"]][["inference"]]
      expect_identical(state[["storage"]], if (mode == "all") "rows" else "cells")
      expect_true(state[["n_cells"]] < n)
      expect_identical(length(state[["resources"]][["views"]]), if (mode == "all") 0L else 1L)
      close(fit)
      expect_identical(DBI::dbListTables(conn), character())
      DBI::dbDisconnect(conn, shutdown = TRUE)
    }
  }
})

local({
  set.seed(1029)
  d = data.frame(x = sample(0:2, 121, TRUE), g = c(rep(1:4, 30), 5),
    h = c(rep(1:3, 40), 4), w = stats::runif(121, 0.3, 1.2),
    c1 = sample(1:9, 121, TRUE), y = stats::rnorm(121))
  conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  fit = dbreg(y ~ x | g + h, conn = conn, data = d, weights = ~w,
    fixef_rm = "singletons", strategy = "compress", retain = "inference")
  reference = fixest::feols(y ~ x | g + h, data = d, weights = ~w,
    fixef.rm = "singletons", ssc = fixest::ssc(K.fixef = "full"))
  expect_equal(stats::nobs(fit), stats::nobs(reference))
  expect_equal(stats::predict(fit), as.numeric(stats::predict(reference)), tolerance = 1e-7)
  expect_equal(as.numeric(stats::vcov(fit, vcov = ~c1)),
    as.numeric(stats::vcov(reference, vcov = ~c1)), tolerance = 1e-7)
  close(fit)
  expect_identical(DBI::dbListTables(conn), character())
  DBI::dbDisconnect(conn, shutdown = TRUE)
})
