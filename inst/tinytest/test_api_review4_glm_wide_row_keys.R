if (requireNamespace("bit64", quietly = TRUE)) local({
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(id = 1:12, x = 1:12, y = 3 + 2*(1:12) + sin(1:12)/10)
  DBI::dbWriteTable(conn, "source", d)
  DBI::dbExecute(conn, paste("CREATE TEMP TABLE keyed AS SELECT *,",
    "CAST('9007199254740992' AS DECIMAL(20,0))+id-1 AS decimal_id,",
    "CAST('9223372036854775808' AS HUGEINT)+id-1 AS huge_id,",
    "CAST('0.10000000000000000000' AS DECIMAL(38,20))+",
    "id*CAST('0.00000000000000000001' AS DECIMAL(38,20)) AS fractional_id FROM source"))
  keys = c("decimal_id", "huge_id", "fractional_id")
  expected = DBI::dbGetQuery(conn, paste("SELECT", paste(paste0("CAST(", keys,
    " AS VARCHAR) AS ", keys), collapse = ","), "FROM keyed ORDER BY id"))
  baseline = sort(DBI::dbListTables(conn))
  expected_fitted = unname(fitted(lm(y ~ x, data = d)))
  check = function(values) {
    for (name in keys) {
      expect_identical(as.character(values[[name]]), expected[[name]])
      expect_identical(length(unique(values[[name]])), nrow(d))
    }
    expect_true(inherits(values$decimal_id, "integer64"))
    expect_true(is.character(values$huge_id))
    expect_true(is.character(values$fractional_id))
  }
  for (maker in list(dbglm)) {
    model = maker(y ~ x, data = d, strategy = "moments")
    result = predict(model, conn = conn, table = "keyed", row_id = keys)
    check(as.data.frame(result))
    expect_equal(as.data.frame(result)$prediction, expected_fitted, tolerance = 1e-10)
    preview = head(result, 3L)
    for (name in keys) expect_identical(as.character(preview[[name]]), expected[[name]][1:3])
    close(result)
    model = maker(y ~ x, conn = conn, table = "keyed", row_id = keys,
      retain = "inference", strategy = "moments")
    saved = list(fitted(model), residuals(model), case.names(model), predict(model, sample = "original"))
    for (rows in saved) check(as.data.frame(rows))
    expect_equal(as.data.frame(saved[[1L]])$.fitted, expected_fitted, tolerance = 1e-10)
    check(generics::augment(model))
    close(model)
    check(as.data.frame(saved[[1L]]))
    for (rows in saved) close(rows)
    expect_identical(sort(DBI::dbListTables(conn)), baseline)
  }
})
