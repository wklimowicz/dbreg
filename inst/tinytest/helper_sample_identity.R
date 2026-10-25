# Shared local/native checks: SQL tuple keys retain their original types and
# names, independently of model expansion, filtering and compression.
dbreg_test_sample_identity = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(id = seq_len(240), x = rep(c(-2, -1, 1, 2), 60),
    g = rep(1:8, each = 30), w = rep(c(.3, .7, 1), 80),
    cl = rep(1:12, each = 20), pair = rep(c("a|b", ""), 120))
  d[["y"]] = 2 + .6 * d[["x"]] + .1 * d[["g"]] + sin(d[["id"]]) / 10
  d[["binary"]] = rep(c(0, 1, 1, 0), 60)
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("sample_identity_"))
  retained = basename(tempfile("sample_retained_"))
  bad = basename(tempfile("sample_bad_"))
  snapshot_resources = NULL
  on.exit({
    dbreg:::release_fit_resources(snapshot_resources)
    dbglm_validation_drop(conn, retained)
    dbglm_validation_drop(conn, bad)
    dbglm_validation_drop(conn, table)
    dbglm_validation_drop(conn, source)
  })
  quote_id = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  keys = c("row key", "dbreg_key_count__")
  fields = quote_id(names(d))
  for (name in c("x", "y", "binary", "cl")) {
    id = switch(name, x = 2, y = 3, binary = 3, cl = 5)
    fields[names(d) == name] = paste0("CASE WHEN id=", id,
      " THEN CAST(NULL AS DOUBLE) ELSE ", quote_id(name), " END AS ", quote_id(name))
  }
  fields[names(d) == "w"] = "CASE WHEN id=4 THEN 0.0 ELSE w END AS w"
  fields = c(fields, paste0("CAST(9007199254740992 AS BIGINT)+CAST(CEIL(id/2.0) AS BIGINT) AS ", quote_id(keys[[1]])),
    paste("pair AS", quote_id(keys[[2]])))
  input_sql = paste("SELECT", paste(fields, collapse = ","), "FROM", source)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, input_sql, c(names(d), keys)))
  inventory = sort(c(before, source, table))
  keep = !d[["id"]] %in% 2:5
  key_values = bit64::as.integer64("9007199254740992") + bit64::as.integer64(ceiling(d[["id"]] / 2))
  expected = data.frame(key = as.character(key_values[keep]), pair = d[["pair"]][keep])
  expected = expected[order(expected[["key"]], expected[["pair"]]), , drop = FALSE]
  rownames(expected) = NULL
  read_keys = function() DBI::dbGetQuery(conn, paste("SELECT", paste(quote_id(keys), collapse = ","),
    "FROM", retained, "ORDER BY", paste(quote_id(keys), collapse = ",")))
  for (strategy in c("moments", "demean", "compress", "mundlak")) {
    fml = if (strategy == "moments") y ~ x else y ~ x | g
    fit = dbreg(fml, conn = conn, table = table, strategy = strategy,
      weights = "w", vcov = ~cl, row_id = keys, sample_table = retained)
    reference = fixest::feols(fml, d[keep, ], weights = ~w, vcov = ~cl,
      fixef.rm = "none", notes = FALSE)
    expect_equal(coef(fit)[["x"]], coef(reference)[["x"]], tolerance = 1e-7, info = strategy)
    expect_equal(nobs(fit), sum(keep), info = strategy)
    expect_identical(fit[["row_id"]], keys)
    expect_identical(fit[["sample_table"]], retained)
    found = read_keys()
    expect_identical(names(found), keys)
    expect_true(inherits(found[[1]], "integer64"))
    comparable = data.frame(key = as.character(found[[1]]), pair = found[[2]])
    expect_identical(comparable, expected, info = strategy)
    expect_identical(dbglm_validation_tables(conn), sort(c(inventory, retained)))
    # A second request cannot replace a caller-owned export.
    expect_error(dbreg(fml, conn = conn, table = table, row_id = keys,
      sample_table = retained), "already exists")
    expect_identical(read_keys(), found)
    dbglm_validation_drop(conn, retained)
  }
  # The same key validation/export layer serves GLM, including unsafe key names
  # and exact BIGINT values. No key is converted through double precision.
  fit = dbglm(binary ~ x, conn = conn, table = table, weights = "w", vcov = ~cl,
    family = binomial("logit"), row_id = keys, sample_table = retained)
  expect_equal(nobs(fit), sum(keep))
  found = read_keys()
  expect_identical(data.frame(key = as.character(found[[1]]), pair = found[[2]]), expected)
  dbglm_validation_drop(conn, retained)
  expect_error(dbglm(binary ~ x, conn = conn, table = table, row_id = keys[[1]],
    sample_table = retained), "unique")
  expect_identical(dbglm_validation_tables(conn), inventory)
  expect_error(dbreg(y ~ x, conn = conn, table = table, row_id = keys[[1]],
    sample_table = retained), "unique")
  expect_identical(dbglm_validation_tables(conn), inventory)
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", table, "LIMIT 0"))
  fields = quote_id(names(schema))
  fields[fields == quote_id(keys[[1]])] = paste0("CASE WHEN id=2 THEN CAST(NULL AS BIGINT) ELSE ",
    quote_id(keys[[1]]), " END AS ", quote_id(keys[[1]]))
  sql = paste("SELECT", paste(fields, collapse = ","), "FROM", table)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, bad, sql, names(schema)))
  # Invalid keys are rejected even on a row otherwise excluded from fitting.
  expect_error(dbreg(y ~ x, conn = conn, table = bad, weights = "w",
    row_id = keys, sample_table = retained), "nonmissing")
  expect_identical(dbglm_validation_tables(conn), sort(c(inventory, bad)))
  dbglm_validation_drop(conn, bad)
  # A validation error AFTER snapshot creation must release that snapshot too.
  fields = quote_id(names(schema))
  fields[names(schema) == "w"] = "CASE WHEN id=6 THEN -0.5 ELSE w END AS w"
  sql = paste("SELECT", paste(fields, collapse = ","), "FROM", table)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, bad, sql, names(schema)))
  expect_error(dbreg(y ~ x, conn = conn, table = bad, weights = "w", row_id = keys,
    sample_table = retained), "negative")
  expect_identical(dbglm_validation_tables(conn), sort(c(inventory, bad)))
  dbglm_validation_drop(conn, bad)
  expect_error(dbreg(y ~ x, conn = conn, table = table, sample_table = retained), "row_id")
  expect_error(dbreg(y ~ x, conn = conn, table = table, row_id = keys,
    sample_table = retained, sql_only = TRUE), "sql_only")
  expect_error(dbreg(y ~ x, conn = conn, table = table, row_id = keys,
    sample_table = retained, data_only = TRUE), "data_only")
  expect_error(dbreg(y ~ x, conn = conn, table = table, row_id = c("id", "id")), "source key columns")
  expect_identical(dbglm_validation_tables(conn), inventory)
  # Snapshot behavior is independent of the caller subsequently replacing its
  # source. This uses real SQL relations, not an R copy or a mocked connection.
  snapshot_resources = dbreg:::fit_resources()
  inputs = new.env(parent = emptyenv())
  inputs[["conn"]] = conn
  inputs[["row_id"]] = keys
  inputs[["sample_table"]] = retained
  inputs[["sample_resources"]] = snapshot_resources
  snapshot_from = dbreg:::prepare_linear_sample(inputs, paste("FROM", table), names(schema))
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, paste(input_sql, "WHERE id>10"),
    names(schema), replace = TRUE))
  expect_equal(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n", snapshot_from))[["n"]]), nrow(d))
  expect_equal(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]]), nrow(d) - 10)
  dbreg:::release_fit_resources(snapshot_resources)
  expect_identical(dbglm_validation_tables(conn), inventory)
  expect_true(DBI::dbIsValid(conn))
}
