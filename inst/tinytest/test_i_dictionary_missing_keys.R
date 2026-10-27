local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  check = function(column, values, expected_levels, reference, expected_name, expected) {
    d = stats::setNames(data.frame(values), column)
    DBI::dbWriteTable(conn, "i_missing_keys", d, overwrite = TRUE)
    formula = stats::as.formula(paste0("~i(", column, ", ref = reference)"))
    compiled = dbreg::sql_model_matrix(formula, conn, "i_missing_keys")
    expect_identical(as.character(compiled[["factor_levels"]][[column]]), expected_levels)
    expect_identical(compiled[["coef_names"]], expected_name)
    expect_identical(anyNA(compiled[["i_specs"]][[1L]][["keys"]]), FALSE)
    sql = paste0("SELECT ", paste(compiled[["select_exprs"]], collapse = ","), " FROM i_missing_keys")
    expect_equal(as.numeric(DBI::dbGetQuery(conn, sql)[[1L]]), expected, tolerance = 1e-14)
  }
  expected = c(0, 1, 0, 0, 1)
  check("f", c("a", "b", NA, "a", "b"), c("a", "b"), "a", "f::b", expected)
  check("f", c(1, 2, NA_real_, 1, 2), c("1", "2"), 1, "f::2", expected)
  check("value", c(2L, 10L, NA_integer_, 2L, 10L), c("2", "10"), 2, "value::10", expected)
  # With no explicit reference, native ordering must choose 2 before 10 even
  # though the outward integer transport is text and its output alias is value.
  DBI::dbWriteTable(conn, "i_missing_keys", data.frame(value = c(2L, 10L, 2L, 10L)), overwrite = TRUE)
  automatic = dbreg::sql_model_matrix(~i(value), conn, "i_missing_keys")
  expect_identical(automatic[["factor_levels"]][["value"]], c("2", "10"))
  expect_identical(automatic[["coef_names"]], "value::10")
  automatic_sql = paste0("SELECT ", automatic[["select_exprs"]], " FROM i_missing_keys")
  expect_equal(as.numeric(DBI::dbGetQuery(conn, automatic_sql)[[1L]]), c(0, 1, 0, 1), tolerance = 1e-14)
  # Existing WHERE clauses stay inside the one-column source. Retaining a
  # missing row in the aliased case must not add a missing category level.
  DBI::dbWriteTable(conn, "i_missing_keys", data.frame(value = c(1L, 2L, 10L, NA_integer_, 2L, 10L)), overwrite = TRUE)
  for (from in c("FROM i_missing_keys WHERE value >= 2",
                 "FROM i_missing_keys AS keys WHERE keys.value >= 2 OR keys.value IS NULL")) {
    filtered = dbreg::sql_model_matrix(~i(value), conn, from)
    expect_identical(filtered[["factor_levels"]][["value"]], c("2", "10"))
    expect_identical(filtered[["coef_names"]], "value::10")
    sql = paste("SELECT", filtered[["select_exprs"]], from)
    expected = if (grepl("IS NULL", from, fixed = TRUE)) c(0, 1, 0, 0, 1) else c(0, 1, 0, 1)
    expect_equal(as.numeric(DBI::dbGetQuery(conn, sql)[[1L]]), expected, tolerance = 1e-14)
  }
  DBI::dbExecute(conn, "CREATE TYPE i_missing_enum AS ENUM ('zeta','alpha')")
  DBI::dbExecute(conn, paste0("CREATE TABLE i_missing_enum_keys AS SELECT CAST(f AS i_missing_enum) AS f FROM ",
    "(VALUES ('zeta'),('alpha'),(NULL),('zeta'),('alpha')) AS source(f)"))
  enum = dbreg::sql_model_matrix(~i(f, ref = "zeta"), conn, "i_missing_enum_keys")
  expect_identical(as.character(enum[["factor_levels"]][["f"]]), c("zeta", "alpha"))
  expect_identical(enum[["coef_names"]], "f::alpha")
  expect_identical(anyNA(enum[["i_specs"]][[1L]][["keys"]]), FALSE)
  enum_sql = paste0("SELECT ", enum[["select_exprs"]], " FROM i_missing_enum_keys")
  expect_equal(as.numeric(DBI::dbGetQuery(conn, enum_sql)[[1L]]), expected, tolerance = 1e-14)
  DBI::dbExecute(conn, paste0("CREATE TABLE i_missing_big_keys AS SELECT CAST(f AS BIGINT) AS f FROM ",
    "(VALUES ('9007199254740992'),('9007199254740993'),(NULL),('9007199254740992'),('9007199254740993')) AS source(f)"))
  big = dbreg::sql_model_matrix(~i(f, ref = "9007199254740992"), conn, "i_missing_big_keys")
  expect_identical(big[["factor_levels"]][["f"]], c("9007199254740992", "9007199254740993"))
  expect_identical(big[["coef_names"]], "f::9007199254740993")
  expect_identical(anyNA(big[["i_specs"]][[1L]][["keys"]]), FALSE)
  big_sql = paste0("SELECT ", big[["select_exprs"]], " FROM i_missing_big_keys")
  expect_equal(as.numeric(DBI::dbGetQuery(conn, big_sql)[[1L]]), expected, tolerance = 1e-14)
  for (values in list(rep(NA_character_, 3), rep(NA_real_, 3), character(), numeric())) {
    DBI::dbWriteTable(conn, "i_missing_keys", data.frame(f = values), overwrite = TRUE)
    expect_error(dbreg::sql_model_matrix(~i(f, ref = "a"), conn, "i_missing_keys"))
  }
  # Nonmissing numeric/character source and frozen selector controls.
  for (values in list(c("a", "b", "c", "a"), c(1.1, 1.2, 1.3, 1.1))) {
    DBI::dbWriteTable(conn, "i_missing_keys", data.frame(f = values), overwrite = TRUE)
    reference = values[[1L]]
    compiled = dbreg::sql_model_matrix(~i(f, ref = reference), conn, "i_missing_keys")
    expect_identical(as.character(compiled[["factor_levels"]][["f"]]), as.character(unique(values)))
    query = paste0("SELECT ", paste(compiled[["select_exprs"]], collapse = ","), " FROM i_missing_keys")
    expected_matrix = cbind(as.numeric(values == values[[2L]]), as.numeric(values == values[[3L]]))
    expect_equal(as.numeric(as.matrix(DBI::dbGetQuery(conn, query))), as.numeric(expected_matrix), tolerance = 1e-14)
  }
})
