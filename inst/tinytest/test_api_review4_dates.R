local({
  d = data.frame(g = rep(as.Date(c("2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04")), each = 20),
    x = rep(rep(-2:2, each = 4), 4))
  d$y = rep(c(1, 10, 100, 1000), each = 20) + .7*d$x + rep(c(-.2, .2), 40)
  nd = data.frame(id = 1:6, g = c("2020-01-01", "bad-date", "2020-01-03", "2020-1-1", "", NA), x = .3)
  expected = c(1.21, NA_real_, 100.21, NA_real_, NA_real_, NA_real_)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "date_points", nd)
  for (maker in list(dbreg)) for (strategy in c("moments", "compress")) {
    for (formula in list(y ~ i(g) + x, y ~ factor(g) + x)) {
      fit = maker(formula, data = d, strategy = strategy)
      expect_equal(unname(predict(fit, nd)), expected, tolerance = 1e-8)
      result = predict(fit, nd, row_id = "id")
      expect_equal(as.data.frame(result)$prediction, expected, tolerance = 1e-8)
      close(result)
      result = predict(fit, conn = conn, table = "date_points", row_id = "id")
      expect_equal(as.data.frame(result)$prediction, expected, tolerance = 1e-8)
      close(result)
      unseen = data.frame(id = 1L, g = "2021-01-01", x = .3)
      if (grepl("factor", deparse(formula), fixed = TRUE)) {
        expect_error(predict(fit, unseen, row_id = "id"), "unseen factor levels")
      } else {
        result = predict(fit, unseen, row_id = "id")
        expect_true(is.na(as.data.frame(result)$prediction))
        close(result)
      }
    }
  }
  old = options(dbreg.prediction_max_entries = getOption("dbreg.prediction_max_entries"))
  on.exit(options(old), add = TRUE)
  for (limit in 2e6) for (maker in list(dbreg)) {
    options(dbreg.prediction_max_entries = limit)
    fit = maker(y ~ x | g, data = d, retain = "prediction",
      strategy = if (identical(maker, dbreg)) "demean" else "moments")
    expect_identical(identical(fit$fitted_state$storage, "database"), limit == 4)
    expect_equal(unname(predict(fit, nd)), expected, tolerance = 1e-8)
    result = predict(fit, nd, row_id = "id")
    expect_equal(as.data.frame(result)$prediction, expected, tolerance = 1e-8)
    close(result)
    close(fit)
  }
})
