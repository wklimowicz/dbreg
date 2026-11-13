if (requireNamespace("duckdb", quietly = TRUE)) {
  set.seed(853)
  d = data.frame(x = rnorm(45), y = rpois(45, 2), extra = letters[rep(1:3, each = 15)])
  d[["date col"]] = as.Date("2024-01-01") + seq_len(nrow(d))
  d$x[4] = NA_real_
  keep = which(!is.na(d$x))
  fits1 = list(dbreg(y ~ x, data = d, retain = "inference"))
  if (exists("dbglm", mode = "function")) {
    fits1 = c(fits1, list(dbglm(y ~ x, data = d, family = poisson(), retain = "inference")))
  }
  for (fit in fits1) {
    result = generics::augment(fit)
    expect_equal(lapply(result[names(d)], identity), lapply(d[keep, , drop = FALSE], identity))
    expect_equal(result$.fitted, fitted(fit))
    expect_equal(result$.resid, residuals(fit, type = "response"))
    supplied = generics::augment(fit, data = d)
    expect_identical(supplied[names(d)], d[keep, , drop = FALSE])
    expect_error(generics::augment(fit, data = d[-1L, ]), "original R input")
    close(fit)
    expect_error(generics::augment(fit), "no live database session")
  }

  if (requireNamespace("bit64", quietly = TRUE)) {
    conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric"))
    DBI::dbWriteTable(conn, "original_input", d)
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE keyed_input AS SELECT *,",
      "CAST('9007199254740992' AS BIGINT) + ROW_NUMBER() OVER () AS id FROM original_input"))
    fits2 = list(dbreg(y ~ x, conn = conn, table = "keyed_input", row_id = "id", retain = "inference"))
    if (exists("dbglm", mode = "function")) {
      fits2 = c(fits2, list(dbglm(y ~ x, conn = conn, table = "keyed_input", row_id = "id",
        family = poisson(), retain = "inference")))
    }
    for (fit in fits2) {
      DBI::dbExecute(conn, "UPDATE keyed_input SET extra = 'changed'")
      result = generics::augment(fit)
      expect_identical(result$extra, d$extra[keep])
      expected = bit64::as.integer64("9007199254740992") + keep
      expect_identical(as.character(result$id), as.character(expected))
      expect_error(generics::augment(fit, data = d), "original R input")
      close(fit)
    }
    expect_identical(sort(DBI::dbListTables(conn)), c("keyed_input", "original_input"))
    DBI::dbDisconnect(conn)
  }
}
