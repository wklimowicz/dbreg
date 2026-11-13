local({
  old = options(dbreg.prediction_max_entries = getOption("dbreg.prediction_max_entries"))
  on.exit(options(old))
  d = expand.grid(g = c("01", "1", "1.0", "1.25", "TRUE", "FALSE", "true", "false"),
    h = c("0", "1"), x = -2:2, rep = 1:4, stringsAsFactors = FALSE)
  d$y = 10*match(d$g, unique(d$g)) + 2*as.numeric(d$h) + .7*d$x + ifelse(d$rep %% 2, -.2, .2)
  for (storage in "portable") {
    options(dbreg.prediction_max_entries = if (storage == "portable") 2e6 else 4)
    for (maker in list(dbreg)) for (formula in list(y ~ x | g, y ~ x | g^h)) {
      fit = maker(formula, data = d, retain = "prediction",
        strategy = if (identical(maker, dbreg)) "demean" else "moments")
      expect_identical(identical(fit$fitted_state$storage, "database"), storage == "database")
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

  # Matrix-produced logical recipes have no explicit key_type. Their SQL
  # training coordinates must retain R's uppercase labels too.
  d = data.frame(g = rep(1:4, 20), flag = rep(c(FALSE, TRUE), each = 40), x = sin(1:80))
  d$y = d$g/3 + 2*d$flag + .7*d$x + cos(1:80)/10
  options(dbreg.prediction_max_entries = 2e6)
  fit = dbreg(y ~ flag + x | g, data = d, strategy = "compress", retain = "prediction")
  oracle = lm(y ~ flag + x + factor(g), data = d)
  expect_true(!is.null(fit$fitted_state))
  expect_equal(predict(fit, d), unname(predict(oracle, d)), tolerance = 1e-8)
  result = predict(fit, transform(d, id = seq_len(nrow(d))), row_id = "id")
  expect_equal(as.data.frame(result)$prediction, unname(predict(oracle, d)), tolerance = 1e-8)
  close(result)
  close(fit)
})
