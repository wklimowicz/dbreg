local({
  # Distinct textual spellings must not let SQL's display format select a
  # different category from the one selected by R's as.character().
  values = c(1, 1.25, 1e-7, 1e20, -1e-7, 1.2345678901234567)
  labels = unique(c("01", as.character(values), "1.0", "0.0000001",
    "100000000000000000000", "1.2345678901234567", "FALSE", "TRUE", "false", "true"))
  d = data.frame(g = rep(labels, each = 20), x = rep(rep(-2:2, each = 4), length(labels)))
  d$y = 10 * rep(seq_along(labels), each = 20) + .7*d$x + rep(c(-.2, .2), nrow(d)/2)
  collect = function(fit, data, ...) {
    result = predict(fit, data, ...)
    on.exit(close(result))
    as.data.frame(result)
  }
  inputs = list(values, c(FALSE, TRUE), c("01", "1", "1.0"))
  for (maker in list(dbreg)) for (strategy in c("moments", "compress")) {
    for (formula in list(y ~ g + x, y ~ i(g) + x, y ~ g*x)) {
      fit = maker(formula, data = d, strategy = strategy)
      for (keys in inputs) {
        nd = data.frame(id = seq_along(keys), g = keys, x = .3)
        expected = 10 * match(as.character(keys), labels) + .21
        expect_equal(unname(predict(fit, nd)), expected, tolerance = 1e-8)
        expect_equal(collect(fit, nd, row_id = "id")$prediction, expected, tolerance = 1e-8)
      }
      # A normalized source may also be a row key. Preserve its native type.
      nd = data.frame(g = c(1.25, 1), x = .3)
      out = collect(fit, nd, row_id = "g")
      expect_identical(out$g, c(1, 1.25))
      expect_equal(out$prediction, 10 * match(as.character(out$g), labels) + .21, tolerance = 1e-8)
    }
  }

  conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "floating_points", data.frame(id = 1:2, g = c(1, 1.25), x = .3))
  DBI::dbWriteTable(conn, "logical_points", data.frame(id = 1:2, g = c(FALSE, TRUE), x = .3))
  for (maker in list(dbreg)) {
    fit = maker(y ~ i(g) + x, data = d, strategy = "moments")
    before = sort(DBI::dbListTables(conn))
    expect_error(predict(fit, conn = conn, table = "floating_points", row_id = "id"),
      "Character prediction keys require")
    expect_identical(sort(DBI::dbListTables(conn)), before)
    result = predict(fit, conn = conn, table = "logical_points", row_id = "id")
    expect_equal(as.data.frame(result)$prediction, 10 * match(c("FALSE", "TRUE"), labels) + .21, tolerance = 1e-8)
    close(result)
  }

  # Exercise the native R VALUES-upload branch with actual local SQL. This
  # controlled branch test does not certify a native Databricks connection.
  native = dbreg:::predict_keyed_sql
  env = new.env(parent = environment(native))
  env$detect_backend = function(conn) {
    info = dbreg:::detect_backend(conn)
    info$is_databricks = TRUE
    info
  }
  env$validate_native_temp_query = function(...) invisible(NULL)
  environment(native) = env
  fit = dbreg(y ~ i(g) + x, data = d, strategy = "moments")
  nd = data.frame(g = c(1.25, 1), x = .3)
  result = native(fit, nd, conn, NULL, NULL, "g", NULL)
  out = as.data.frame(result)
  expect_identical(out$g, c(1, 1.25))
  expect_equal(out$prediction, 10 * match(as.character(out$g), labels) + .21, tolerance = 1e-8)
  close(result)

  # SQL typeof(), rather than the driver's numeric BIGINT collection mode,
  # determines whether conversion to a character key preserves exactness.
  keys = c("-9223372036854775808", "9007199254740993", "9223372036854775807")
  d = data.frame(g = rep(keys, each = 20), y = rep(c(1, 10, 100), each = 20) + rep((1:20)/100, 3))
  DBI::dbExecute(conn, paste0("CREATE TEMP TABLE integer_points AS SELECT * FROM (VALUES ",
    paste0("(", seq_along(keys), ",CAST('", keys, "' AS BIGINT))", collapse = ","), ") p(id,g)"))
  for (maker in list(dbreg)) {
    fit = maker(y ~ i(g), data = d, strategy = "moments")
    result = predict(fit, conn = conn, table = "integer_points", row_id = "id")
    expect_equal(as.data.frame(result)$prediction, c(1.105, 10.105, 100.105), tolerance = 1e-8)
    close(result)
  }
})
