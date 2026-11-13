# Numeric-looking character labels remain distinct under SQL coercion. SQL
# predictions must agree with the R interpretation of the learned key type.
local({
  old = options(dbreg.prediction_max_entries = getOption("dbreg.prediction_max_entries"))
  on.exit(options(old), add = TRUE)
  d = data.frame(id = 1:80, g = rep(c("0", "01", "1", "2"), each = 20))
  d[["y"]] = rep(c(0, 10, 100, 1000), each = 20) + rep((1:20)/100, 4)
  for (fml in list(y ~ g, y ~ i(g))) for (strategy in c("moments", "compress")) {
    fit = dbreg(fml, data = d, strategy = strategy)
    for (keys in list(c(1L, 2L, 0L), c("01", "1", "2", "0"))) {
      newdata = data.frame(id = seq_along(keys), g = keys)
      prediction = predict(fit, newdata, row_id = "id")
      got = as.data.frame(prediction)[[prediction[["prediction_column"]]]]
      expect_equal(got, unname(predict(fit, newdata)), tolerance = 1e-10)
      if (identical(keys[[1L]], 1L)) expect_equal(got[[1L]], 100.105, tolerance = 1e-10)
      close(prediction)
    }
  }
  # Exercise portable FE dictionaries and composite keys.
  # Session-backed dictionaries have additional coverage in their own commit.
  for (dictionary_limit in 1e6) {
    saved = options(dbreg.prediction_max_entries = dictionary_limit)
    d[["x"]] = sin(d[["id"]])
    d[["h"]] = rep(c("01", "1"), 40)
    d[["y"]] = d[["y"]] + .5*d[["x"]]
    for (fml in list(y ~ x | g, y ~ x | g^h)) {
      fit = dbreg(fml, data = d, retain = "prediction", strategy = "demean")
      nd = data.frame(id = 1:2, g = c(1L, 2L), h = c(1L, 1L), x = c(.1, .2))
      p = predict(fit, nd, row_id = "id")
      expect_equal(as.data.frame(p)[[p[["prediction_column"]]]], unname(predict(fit, nd)), tolerance = 1e-8)
      close(p)
      close(fit)
    }
    options(saved)
  }
})
