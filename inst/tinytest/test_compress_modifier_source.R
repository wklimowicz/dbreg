# The same raw numeric values can supply FE, structural, weight and offset roles.
local({
  set.seed(78524)
  n = 120L
  d = data.frame(g = rep(c(10L, 20L, 40L), length.out = n), x = rnorm(n))
  d[["y"]] = d[["g"]] + .4 * d[["x"]] + sin(d[["g"]]) + rnorm(n, sd = .2)
  X = model.matrix(~x + factor(g), d)
  ref = lm.wfit(X, d[["y"]] - d[["g"]], d[["g"]])
  bread = solve(crossprod(X, d[["g"]] * X))
  score = X * (d[["g"]] * ref[["residuals"]])
  V = bread %*% crossprod(score) %*% bread * n / (n - ncol(X))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  before = sort(DBI::dbListTables(conn))
  model = dbreg(y ~ x + g | g, conn = conn, data = d, strategy = "compress",
    weights = "g", offset = "g", vcov = "hc1", retain = "none")
  expect_identical(names(coef(model)), "x")
  expect_equal(as.numeric(coef(model)), as.numeric(ref[["coefficients"]][2L]), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(model)), as.numeric(V[2L, 2L]), tolerance = 1e-8)
  expect_equal(attr(vcov(model), "df.K"), ncol(X))
  expect_true(is.na(coef(model, collin = TRUE)[["g"]]))
  expected = as.numeric(ref[["fitted.values"]] + d[["g"]])
  expect_equal(predict(model, d), expected, tolerance = 1e-8)
  minimal = d[c("x", "g")]
  expect_equal(predict(model, minimal), expected, tolerance = 1e-8)
  minimal[["g"]][2L] = NA_integer_
  minimal[["x"]][5L] = NA_real_
  expected[c(2L, 5L)] = NA_real_
  actual = predict(model, minimal)
  expect_identical(is.na(actual), is.na(expected))
  expect_equal(actual, expected, tolerance = 1e-8)
  close(model)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), before)
})
