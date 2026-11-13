local({
  set.seed(11150)
  d = data.frame(key = sprintf("id-%03d", 1:100),
    g = rep(1:4, each = 25), x = rep(1:5, 20), y = rnorm(100), w = runif(100,.2,4))
  d[["x"]][c(2L, 9L)] = NA
  keep = !is.na(d[["x"]])
  ref = glm(y ~ x + factor(g), data = d, weights = w)
  model = dbglm(y ~ x | g, data = d, weights = ~w, retain = "all", strategy = "compress")
  expect_equal(fitted(model), unname(fitted(ref)), tolerance = 1e-10)
  expect_equal(predict(model, type = "link"), unname(fitted(ref)), tolerance = 1e-10)
  expect_identical(case.names(model), which(keep))
  for (type in c("response", "deviance", "pearson", "working")) {
    expect_equal(residuals(model, type), unname(residuals(ref, type = type)), tolerance = 1e-10)
  }
  expect_identical(which(is.na(residuals(model, na.rm = FALSE))), which(!keep))
  owned = model[["resource_state"]][["conn"]]
  close(model)
  close(model)
  expect_false(DBI::dbIsValid(owned))
  expect_error(fitted(model), "no live database session")

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "retained_src", d)
  model = dbglm(y ~ x | g, conn = conn, table = "retained_src", row_id = "key",
    weights = ~w, retain = "inference", strategy = "compress")
  expect_equal(length(DBI::dbListTables(conn)), 3L)
  lazy = fitted(model)
  expect_true(inherits(lazy, "dbreg_glm_rows"))
  cases = case.names(model)
  expect_identical(as.data.frame(cases)[["key"]], d[["key"]][keep])
  close(cases)
  result = as.data.frame(lazy)
  expect_identical(result[["key"]], d[["key"]][keep])
  expect_equal(result[[".fitted"]], unname(fitted(ref)), tolerance = 1e-10)
  residual_result = residuals(model, na.rm = FALSE)
  residual = as.data.frame(residual_result)
  close(residual_result)
  expect_identical(residual[["key"]], d[["key"]])
  expect_identical(which(is.na(residual[[".residual"]])), which(!keep))
  close(model)
  expect_equal(as.data.frame(lazy), result)
  expect_equal(head(lazy, 3L), utils::head(result, 3L))
  close(lazy)
  close(lazy)
  expect_error(as.data.frame(lazy), "no live database session")
  expect_true(DBI::dbIsValid(conn))
  expect_identical(DBI::dbListTables(conn), "retained_src")
  model = dbglm(y ~ x | g, data = d, row_id = "key", retain = "all")
  child = fitted(model)
  own = model[["resource_state"]][["conn"]]
  close(model)
  expect_true(DBI::dbIsValid(own))
  expect_identical(as.data.frame(child)[["key"]], d[["key"]][keep])
  close(child)
  expect_false(DBI::dbIsValid(own))

  if (requireNamespace("bit64", quietly = TRUE)) {
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE exact_keys AS SELECT",
      "CAST('9007199254740992' AS BIGINT) + ROW_NUMBER() OVER () AS id,x,y,w FROM retained_src"))
    model = dbglm(y ~ x, conn = conn, table = "exact_keys", row_id = "id",
      weights = ~w, retain = "all")
    lazy = fitted(model)
    close(model)
    keys = as.data.frame(lazy)[["id"]]
    expect_identical(length(unique(keys)), sum(keep))
    expected = bit64::as.integer64("9007199254740992") + seq_len(nrow(d))
    expect_identical(as.character(keys), as.character(expected[keep]))
    close(lazy)
    DBI::dbExecute(conn, "DROP TABLE exact_keys")
  }

  expect_error(dbglm(y ~ x | g, conn = conn, table = "retained_src", row_id = "g",
    retain = "inference"), "unique and nonmissing")
  expect_identical(DBI::dbListTables(conn), "retained_src")
})
