if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("bit64", quietly = TRUE)) {
  set.seed(98)
  d = expand.grid(g = c("north:1", "south"), h = 1:2, x = -1:1, trend = -1:1, rep = 1:4)
  d$id = seq_len(nrow(d))
  d$h = bit64::as.integer64("9007199254740992") + d$h
  d$o = d$rep / 20
  d$y = rpois(nrow(d), exp(.7 + d$x / 5 + d$trend / 4 + d$o))
  models = list(dbreg(y ~ x | (g^h)[trend], data = d, offset = ~o,
    strategy = "compress", retain = "prediction"),
    dbglm(y ~ x | (g^h)[trend], data = d, offset = ~o, family = poisson(),
      strategy = "compress", retain = "prediction"))
  pred = d[c("id", "g", "h", "x", "trend", "o")]
  pred$h[5L] = bit64::as.integer64("9007199254740999")
  pred$trend[12L] = NA_real_
  for (fit in models) {
    close(fit)
    ref = predict(fit, pred)
    expect_identical(which(is.na(ref)), c(5L, 12L))
    sql = predict(fit, pred, row_id = "id")
    expect_equal(as.data.frame(sql)[[sql$prediction_column]], ref, tolerance = 1e-10)
    close(sql)
    previous = options(dbreg.projection_max_upload_bytes = 1)
    expect_error(predict(fit, pred, row_id = "id"), "upload budget")
    options(previous)
  }
  d = data.frame(g = rep(c("a", "b"), each = 6), h = rep(c("c", "d"), each = 6),
    x = rep(-1:1, 4), y = rnorm(12), id = 1:12)
  d$copy = d$x
  model = dbreg(y ~ x + copy | g + h, data = d, strategy = "demean", retain = "prediction")
  pred = d[1:3, ]
  pred$h[2] = "d"
  pred$copy[3] = 4
  result = predict(model, pred, row_id = "id")
  close(model)
  expect_equal(as.data.frame(result)[[result$prediction_column]], predict(model, pred), tolerance = 1e-10)
  expect_identical(which(is.na(as.data.frame(result)[[result$prediction_column]])), c(2L, 3L))
  close(result)
}
