# Session-backed learned keys obey the same character identity as portable
# dictionaries, even when the prediction input has a numeric SQL type.
local({
  saved = options(dbreg.prediction_max_entries = 4)
  on.exit(options(saved), add = TRUE)
  d = data.frame(id = 1:80, g = rep(c("0", "01", "1", "2"), each = 20),
    h = rep(c("01", "1"), 40), x = sin(1:80))
  d[["y"]] = rep(c(0, 10, 100, 1000), each = 20) + rep((1:20)/100, 4) + .5*d[["x"]]
  for (fml in list(y ~ x | g, y ~ x | g^h)) {
    fit = dbreg(fml, data = d, retain = "prediction", strategy = "demean")
    expect_identical(fit[["fitted_state"]][["storage"]], "database")
    nd = data.frame(id = 1:2, g = c(1L, 2L), h = c(1L, 1L), x = c(.1, .2))
    p = predict(fit, nd, row_id = "id")
    expect_equal(as.data.frame(p)[[p[["prediction_column"]]]], unname(predict(fit, nd)), tolerance = 1e-8)
    close(p)
    close(fit)
  }
})
