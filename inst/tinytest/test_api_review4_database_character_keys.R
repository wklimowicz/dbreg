local({
  old = options(dbreg.prediction_max_entries = getOption("dbreg.prediction_max_entries"))
  on.exit(options(old))
  d = expand.grid(g = c("01", "1", "1.0", "1.25", "TRUE", "FALSE", "true", "false"),
    h = c("0", "1"), x = -2:2, rep = 1:4, stringsAsFactors = FALSE)
  d$y = 10*match(d$g, unique(d$g)) + 2*as.numeric(d$h) + .7*d$x + ifelse(d$rep %% 2, -.2, .2)
  for (storage in "database") {
    options(dbreg.prediction_max_entries = if (storage == "portable") 2e6 else 4)
    for (maker in list(dbreg)) for (formula in list(y ~ x | g, y ~ x | g^h)) {
      fit = maker(formula, data = d, retain = "prediction",
        strategy = if (identical(maker, dbreg)) "demean" else "moments")
      expect_identical(dbreg:::is_database_effect_state(fit$fitted_state), storage == "database")
      for (keys in list(c(1.25, 1), c(TRUE, FALSE))) {
        nd = data.frame(id = 1:2, g = keys, h = c(1L, 0L), x = c(.3, -.2))
        # The one-family model absorbs the common h contribution at its mean.
        expected = 10*match(as.character(keys), unique(d$g)) +
          if (grepl("^", deparse(formula), fixed = TRUE)) 2*nd$h + .7*nd$x else 1 + .7*nd$x
        expect_equal(unname(predict(fit, nd)), expected, tolerance = 1e-8)
        result = predict(fit, nd, row_id = "id")
        expect_equal(as.data.frame(result)$prediction, expected, tolerance = 1e-8)
        close(result)
      }
      close(fit)
    }
  }

})
