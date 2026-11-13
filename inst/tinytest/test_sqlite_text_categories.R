if (!requireNamespace("RSQLite", quietly = TRUE)) exit_file("RSQLite not installed")
local({
  conn = DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(conn))
  d = data.frame(x = rep(c(-1, 0, 2, 3), 9), g = rep(c("a", "b", "c"), each = 12))
  d[["y"]] = 1 + .7 * d[["x"]] + match(d[["g"]], c("a", "b", "c")) +
    .2 * sin(seq_len(nrow(d)))
  DBI::dbWriteTable(conn, "text_categories", d)
  initial = sort(DBI::dbListTables(conn))
  reference = stats::lm(y ~ x + factor(g), data = d)
  for (strategy in c("moments", "compress")) {
    for (fml in list(y ~ x + factor(g), y ~ x + i(g))) {
      fit = dbreg(fml, conn = conn, table = "text_categories", strategy = strategy)
      expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-9)
    }
    # A nondefault reference also takes the typed-dictionary route for bare text.
    fit = dbreg(y ~ x + g, conn = conn, table = "text_categories", strategy = strategy,
      factor_refs = list(g = "c"))
    input = d
    input[["g"]] = stats::relevel(factor(input[["g"]]), "c")
    expected = stats::lm(y ~ x + g, data = input)
    expect_equal(unname(coef(fit)), unname(coef(expected)), tolerance = 1e-9)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(expected)), tolerance = 1e-9)
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # SQLite REAL is 64-bit IEEE-754 double, and a NULL in row 1 must not hide the column type.
  dn = data.frame(x = c(0, rep(c(-1, 0, 2, 3), 9)), g = c(NA_real_, rep(c(0.5, 1.5, 2.5), each = 12)))
  dn[["y"]] = 1 + .7 * dn[["x"]] + ifelse(is.na(dn[["g"]]), 0, dn[["g"]]) + .2 * sin(seq_len(nrow(dn)))
  DBI::dbWriteTable(conn, "real_categories", dn)
  ref_n = stats::lm(y ~ x + factor(g), data = dn[-1L, ])
  for (strategy in c("moments", "compress")) {
    for (fml in list(y ~ x + factor(g), y ~ x + i(g))) {
      fit = dbreg(fml, conn = conn, table = "real_categories", strategy = strategy)
      expect_equal(unname(coef(fit)), unname(coef(ref_n)), tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref_n)), tolerance = 1e-9)
    }
  }
  mm = dbreg:::sql_model_matrix(~ x + i(g), conn = conn, table = "real_categories")
  expect_equal(length(mm[["col_names"]]), 3L)
})

## Logical and Date fixed effects: SQLite stores them as numbers, so R logical/Date
## prediction data must error rather than silently miss every level; DuckDB keeps the types.
if (requireNamespace("fixest", quietly = TRUE)) local({
  set.seed(1)
  n = 600
  d = data.frame(x = rnorm(n), union = runif(n) > .6, firm = sample(20, n, TRUE),
    month = as.Date("2015-01-01") + 30 * sample(0:11, n, TRUE))
  d$y = 1 + .5 * d$x + .3 * d$union + rnorm(20)[d$firm] + rnorm(n)
  sq = DBI::dbConnect(RSQLite::SQLite())
  DBI::dbWriteTable(sq, "d", d)
  fit = dbreg(y ~ x | union + firm, conn = sq, table = "d", strategy = "compress")
  expect_error(predict(fit, d), "union")
  stored = transform(d, union = as.integer(union))
  expect_equal(unname(predict(fit, stored)), unname(predict(fixest::feols(y ~ x | union + firm, d), d)),
    tolerance = 1e-10)
  fit = dbreg(y ~ x | month, conn = sq, table = "d", strategy = "compress")
  expect_error(predict(fit, d), "month")
  DBI::dbDisconnect(sq)
  for (f in list(y ~ x | union + firm, y ~ x | month)) {
    fit = dbreg(f, data = d, strategy = "compress")
    expect_equal(unname(predict(fit, d)), unname(predict(fixest::feols(f, d), d)), tolerance = 1e-10)
  }
})
