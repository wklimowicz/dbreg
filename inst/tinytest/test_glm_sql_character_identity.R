# The GLM structural prediction compiler must preserve learned character keys.
local({
  d = data.frame(id = 1:80, g = rep(c("0", "01", "1", "2"), each = 20))
  d[["y"]] = rep(c(0, 10, 100, 1000), each = 20) + rep((1:20)/100, 4)
  for (strategy in c("moments", "compress")) {
    fit = dbglm(y ~ g, data = d, strategy = strategy)
    for (keys in list(c(1L, 2L, 0L), c("01", "1", "2", "0"))) {
      nd = data.frame(id = seq_along(keys), g = keys)
      p = predict(fit, nd, row_id = "id")
      expect_equal(as.data.frame(p)[[p[["prediction_column"]]]], unname(predict(fit, nd)), tolerance = 1e-10)
      close(p)
    }
  }
})
