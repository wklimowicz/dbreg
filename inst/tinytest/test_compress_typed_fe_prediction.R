# Typed i() expansion reads raw keys before replaying a shared categorical FE.
local({
  set.seed(56419)
  n = 120L
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  original_tables = sort(DBI::dbListTables(conn))
  inputs = list(
    list(keys = c(1, 1 + 2 * .Machine[["double.eps"]], 2), decimal = FALSE),
    list(keys = c(FALSE, TRUE), decimal = FALSE),
    list(keys = c(10L, 20L, 40L), decimal = TRUE)
  )
  for (input in inputs) {
    keys = input[["keys"]]
    d = data.frame(x = rnorm(n), g = rep(keys, length.out = n))
    group = match(d[["g"]], keys)
    D = diag(length(keys))[group, , drop = FALSE]
    X = cbind(1, d[["x"]], D[, -1L, drop = FALSE])
    d[["y"]] = .6 * d[["x"]] + group + rnorm(n, sd = .2)
    reference = lm.fit(X, d[["y"]])
    model = dbreg(y ~ x + i(g) | g, conn = conn, data = d,
      strategy = "compress", vcov = "hc1", retain = "none")
    expect_equal(as.numeric(coef(model)), as.numeric(reference[["coefficients"]][2L]), tolerance = 1e-8)
    expect_equal(attr(vcov(model), "df.K"), ncol(X))
    expect_equal(predict(model, d), as.numeric(reference[["fitted.values"]]), tolerance = 1e-8)
    minimal = d[c("x", "g")]
    if (input[["decimal"]]) minimal[["g"]] = as.character(minimal[["g"]])
    expect_equal(predict(model, minimal), as.numeric(reference[["fitted.values"]]), tolerance = 1e-8)
    expect_equal(predict(model, minimal[c(9L, 2L, 7L), , drop = FALSE]),
      as.numeric(reference[["fitted.values"]][c(9L, 2L, 7L)]), tolerance = 1e-8)
    unseen = minimal[1L, , drop = FALSE]
    unseen[["g"]] = if (is.logical(keys)) "not logical" else if (input[["decimal"]]) "99" else 1 + .Machine[["double.eps"]]
    expect_error(predict(model, unseen))
    close(model)
  }
  # A categorical typed first also supplies ordinary factor terms. A one-level
  # prediction batch must keep the full saved labels and reference coding.
  d = expand.grid(f = c("a", "b", "c", "d"), h = seq_len(5L), repeat_id = seq_len(6L),
    stringsAsFactors = FALSE)
  d[["x"]] = rnorm(n)
  X = model.matrix(~x + f + factor(h), d)
  d[["y"]] = as.numeric(X %*% c(.2, .7, .3, -.4, .8, seq_len(4L) / 5)) + rnorm(n, sd = .2)
  reference = lm.fit(X, d[["y"]])
  model = dbreg(y ~ x + f + i(f, keep = c("b", "d")) | h, conn = conn, data = d,
    strategy = "compress", retain = "none")
  expect_equal(as.numeric(coef(model)), as.numeric(reference[["coefficients"]][2:5]), tolerance = 1e-8)
  expect_equal(predict(model, d), as.numeric(reference[["fitted.values"]]), tolerance = 1e-8)
  rows = which(d[["f"]] == "d")[c(3L, 1L, 7L)]
  minimal = d[rows, c("x", "f", "h"), drop = FALSE]
  expect_equal(predict(model, minimal), as.numeric(reference[["fitted.values"]][rows]), tolerance = 1e-8)
  close(model)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), original_tables)
})
