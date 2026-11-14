# A constant FE is an intercept, whether it is the sole FE or a redundant one.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(id = 1:120, x = rep(1:3, 40), g = "only",
    h = rep(letters[1:4], each = 30), w = rep(c(.5, 1, 2, 3), 30),
    cl = rep(1:12, 10))
  d[["y"]] = 2 + .7 * d[["x"]] + .3 * match(d[["h"]], letters[1:4]) + sin(d[["id"]])
  for (second in c(FALSE, TRUE)) for (intercept in c(FALSE, TRUE)) {
    fml = if (second) {
      if (intercept) y ~ x | g + h else y ~ 0 + x | g + h
    } else if (intercept) y ~ x | g else y ~ 0 + x | g
    ref = stats::lm(if (second) y ~ x + h else y ~ x, data = d, weights = w)
    X = stats::model.matrix(ref)
    bread = solve(crossprod(X, X * d[["w"]]))
    scores = X * (d[["w"]] * residuals(ref))
    for (strategy in c("auto", "compress")) for (kind in c("iid", "hc1", "cluster")) {
      request = if (kind == "cluster") ~cl else kind
      fit = dbreg(fml, data = d, conn = conn, weights = "w", vcov = request,
        strategy = strategy, retain = "all")
      expect_identical(fit[["strategy"]], "compress")
      expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-9)
      meat = if (kind == "cluster") crossprod(rowsum(scores, d[["cl"]])) else crossprod(scores)
      K = ncol(X)
      expected = if (kind == "iid") vcov(ref) else bread %*% meat %*% bread *
        if (kind == "cluster") 12/11 * (nrow(d)-1)/(nrow(d)-K) else nrow(d)/(nrow(d)-K)
      expect_equal(as.numeric(vcov(fit)), as.numeric(expected["x", "x"]), tolerance = 1e-9)
      expect_equal(fitted(fit), unname(fitted(ref)), tolerance = 1e-9)
      expect_equal(predict(fit, d), unname(fitted(ref)), tolerance = 1e-9)
      expect_equal(df.residual(fit), df.residual(ref))
      unseen = d[1:2, ]
      unseen[["g"]] = c("other", NA_character_)
      expect_true(all(is.na(predict(fit, unseen))))
      close(fit)
    }
  }
  for (fml in list(y ~ x | g, y ~ 0 + x | g, y ~ x | g + h)) {
    fit = dbreg(fml, data = d, conn = conn, strategy = "compress")
    expect_true(all(is.finite(predict(fit))))
    reference = stats::lm(if (length(fit[["fe"]]) > 1L) y ~ x + h else y ~ x, data = d)
    expect_equal(predict(fit, d), unname(fitted(reference)), tolerance = 1e-9)
  }
  # Numeric one-level keys have the same constant-FE interpretation.
  d[["g"]] = 1L
  fit = dbreg(y ~ x | g, data = d, conn = conn, strategy = "compress", retain = "all")
  expect_equal(fitted(fit), unname(fitted(stats::lm(y ~ x, data = d))), tolerance = 1e-9)
  close(fit)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
