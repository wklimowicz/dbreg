# Poisson's fixed dispersion does not require positive residual degrees of freedom.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(y = c(2, 4, 6), g = letters[1:3], w = c(.5, 2, 3))
  for (strategy in c("moments", "compress")) for (weighted in c(FALSE, TRUE)) {
    weights = if (weighted) "w" else NULL
    expected = diag(1 / (d[["y"]] * if (weighted) d[["w"]] else 1))
    for (inference in c(FALSE, TRUE)) {
      fit = dbglm(y ~ 0 + g, data = d, conn = conn, family = poisson(),
        strategy = strategy, weights = weights, ssc = list(K.adj = FALSE),
        control = dbglm_control(compute_inference = inference), retain = "inference")
      expect_true(fit[["converged"]] && fit[["finite_mle"]])
      expect_equal(unname(coef(fit)), log(d[["y"]]), tolerance = 1e-9)
      expect_equal(predict(fit, d, type = "response"), d[["y"]], tolerance = 1e-9)
      expect_equal(deviance(fit), 0, tolerance = 1e-12)
      expect_equal(df.residual(fit), 0)
      if (inference) {
        expect_equal(as.numeric(vcov(fit)), as.numeric(expected), tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit, vcov = "iid")), as.numeric(expected), tolerance = 1e-8)
        expect_true(all(is.finite(confint(fit))))
        expect_error(vcov(fit, vcov = "iid", ssc = list(K.adj = TRUE)), "degrees of freedom")
      } else expect_identical(fit[["inference_result"]][["reason"]], "not_requested")
      close(fit)
    }
    expect_error(dbglm(y ~ 0 + g, data = d, conn = conn, family = poisson(),
      strategy = strategy), "degrees of freedom")
    expect_error(dbglm(y ~ 0 + g, data = d, conn = conn, family = gaussian(),
      strategy = strategy, ssc = list(K.adj = FALSE),
      control = dbglm_control(compute_inference = FALSE)), "degrees of freedom")
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
