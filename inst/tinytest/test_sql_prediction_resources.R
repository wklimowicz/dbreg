local({
  original = options(dbreg.prediction_max_entries = 64)
  on.exit(options(original))
  set.seed(942)
  d = data.frame(id = 1:240, g = rep(1:24, each = 10), x = rnorm(240), y = rnorm(240))
  fit = dbreg(y ~ x | g, data = d, retain = "prediction", strategy = "demean")
  output = predict(fit, d[1:12, ], row_id = "id")
  conn = output[["conn"]]
  compact = dbreg(y ~ x, data = d, strategy = "moments")
  output2 = predict(compact, newdata = d[13:22, ], conn = conn, row_id = "id")
  expected = as.data.frame(output2)
  expect_error(predict(compact, d[1:3, ], conn = conn, row_id = "id", prediction_table = "private_export"), "caller-owned")
  close(fit)
  expect_true(DBI::dbIsValid(conn))
  close(output)
  expect_true(DBI::dbIsValid(conn))
  expect_equal(as.data.frame(output2), expected)
  # Explicitly passing a released public handle must not bypass its lease by
  # replacing it with another live owner from the fitted model.
  borrowed = dbreg(y ~ x | g, conn = output2[["conn"]], data = d, strategy = "demean", retain = "prediction")
  expect_error(predict(borrowed, d[1:3, ], conn = conn, row_id = "id"), "released")
  close(borrowed)
  expect_error(predict(compact, d[1:3, ], conn = conn, row_id = "id"), "released")
  close(output2)
  expect_false(DBI::dbIsValid(conn))
  close(output2)

  # Passing a live public connection to another fit also shares its ownership.
  result = predict(compact, d[1:12, ], row_id = "id")
  shared = result[["conn"]]
  retained = dbreg(y ~ x | g, conn = shared, data = d, strategy = "demean", retain = "prediction")
  expected = predict(retained, d[1:3, ])
  close(result)
  expect_true(DBI::dbIsValid(shared))
  expect_equal(predict(retained, d[1:3, ]), expected)
  close(retained)
  expect_false(DBI::dbIsValid(shared))

  # An abandoned output's finalizer drops its table, while caller sessions live.
  caller = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(caller, shutdown = TRUE), add = TRUE)
  abandoned = predict(compact, d[1:3, ], conn = caller, row_id = "id")
  name = abandoned[["table"]]
  rm(abandoned)
  invisible(gc())
  expect_false(DBI::dbExistsTable(caller, name))
  expect_true(DBI::dbIsValid(caller))
})
