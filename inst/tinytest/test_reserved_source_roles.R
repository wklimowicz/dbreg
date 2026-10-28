# Original SQL names remain quoted at the source boundary; private aliases are
# used downstream. Compare public fits independently of identifier rewriting.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(x = rep(c(-1.5, -.5, .5, 1.5), 30))
  d[["row key"]] = seq_len(nrow(d))
  d[["other field"]] = sin(seq_len(nrow(d)))
  d[["binary"]] = as.numeric(seq_len(nrow(d)) %% 7L %in% c(0L, 2L, 3L))
  d[["response"]] = as.numeric(seq_len(nrow(d)) %% 11L %in% c(0L, 1L, 4L, 6L, 8L))
  d[["y"]] = 1 + .6 * d[["binary"]] + .3 * d[["other field"]] + cos(seq_len(nrow(d))) / 10
  DBI::dbWriteTable(conn, "reserved_source_roles", d)
  source = dbreg:::safe_source(conn, "FROM reserved_source_roles", c("binary", "x", "row key"))
  expect_false(identical(source[["to_safe"]]("binary"), "binary"))
  expect_identical(source[["to_safe"]]("x"), "x")
  expect_true(grepl(as.character(DBI::dbQuoteIdentifier(conn, "x")), source[["from_statement"]], fixed = TRUE))
  expect_identical(nrow(DBI::dbGetQuery(conn, paste("SELECT *", source[["from_statement"]]))), nrow(d))
  shuffled = d[rev(seq_len(nrow(d))), , drop = FALSE]
  for (fml in list(binary ~ x + `other field`, y ~ binary + `other field`)) {
    fit = dbreg(fml, conn = conn, table = "reserved_source_roles", row_id = "row key", strategy = "moments")
    reference = stats::lm(fml, data = d)
    expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-8)
    expect_identical(dim(vcov(fit)), dim(vcov(reference)))
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-8)
    expect_equal(as.numeric(predict(fit, newdata = shuffled)), as.numeric(predict(reference, newdata = shuffled)), tolerance = 1e-8)
    expect_equal(nobs(fit), nrow(d))
  }
  for (fml in list(binary ~ x + `other field`, response ~ binary + `other field`)) {
    fit = dbglm(fml, conn = conn, table = "reserved_source_roles", row_id = "row key",
      family = stats::binomial(), ssc = list(K.adj = FALSE))
    reference = stats::glm(fml, data = d, family = stats::binomial())
    expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-7)
    expect_identical(dim(vcov(fit)), dim(vcov(reference)))
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-7)
    expect_equal(as.numeric(predict(fit, newdata = shuffled, type = "response")), as.numeric(predict(reference, newdata = shuffled, type = "response")), tolerance = 1e-7)
    expect_equal(nobs(fit), nrow(d))
  }
})
