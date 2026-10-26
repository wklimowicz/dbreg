local({
  set.seed(13013)
  d = data.frame(grp = rep(c("a", "b"), 120), grpb = stats::rnorm(240))
  d[["y"]] = 2 + 3 * (d[["grp"]] == "b") - 0.8 * d[["grpb"]] + stats::rnorm(240)
  reference = stats::lm(y ~ grp + grpb, d)
  design = stats::model.matrix(reference)
  inverse = solve(crossprod(design))
  scores = design * as.numeric(stats::residuals(reference))
  covariance = inverse %*% crossprod(scores) %*% inverse * nrow(d) / (nrow(d) - ncol(design))
  intervals = function(V, df) {
    error = stats::qt(0.975, df) * sqrt(diag(V))
    cbind(as.numeric(stats::coef(reference)) - error, as.numeric(stats::coef(reference)) + error)
  }
  expected_iid = intervals(stats::vcov(reference), stats::df.residual(reference))
  expected = intervals(covariance, stats::df.residual(reference))
  fit = dbreg::dbreg(y ~ grp + grpb, data = d, strategy = "moments", retain = "inference")
  on.exit(close(fit), add = TRUE)
  changed = summary(fit, vcov = "hc1")
  expect_equal(as.numeric(stats::confint(fit)), as.numeric(expected_iid), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(changed)), as.numeric(expected), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(fit, vcov = "hc1")), as.numeric(expected), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(changed, parm = c(3, 2, 3))), as.numeric(expected[c(3, 2, 3), ]), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(changed, parm = -1)), as.numeric(expected[-1, ]), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(changed, parm = c(0, 3, 2))), as.numeric(expected[c(0, 3, 2), ]), tolerance = 1e-10)
  expect_equal(as.numeric(stats::confint(changed, parm = "grpb")), as.numeric(expected[2, , drop = FALSE]), tolerance = 1e-10)
  expect_identical(dim(stats::confint(changed, parm = numeric())), c(0L, 2L))
  expect_error(stats::confint(changed, parm = 4), "Unknown coefficient")
  expect_error(stats::confint(changed, parm = "absent"), "Unknown coefficient")

  for (model in list(fit, changed)) {
    table = generics::tidy(model, conf.int = TRUE)
    ci = if (isTRUE(model[["postfit_covariance"]])) expected else expected_iid
    expect_equal(table[["estimate"]], as.numeric(stats::coef(reference)), tolerance = 1e-10)
    expect_equal(table[["conf.low"]], as.numeric(ci[, 1]), tolerance = 1e-10)
    expect_equal(table[["conf.high"]], as.numeric(ci[, 2]), tolerance = 1e-10)
    full_table = generics::tidy(model, collin = TRUE, conf.int = TRUE)
    expect_equal(full_table[["estimate"]], as.numeric(stats::coef(reference)), tolerance = 1e-10)
    expect_equal(full_table[["conf.low"]], as.numeric(ci[, 1]), tolerance = 1e-10)
  }

  attr(covariance, "df.t") = 7
  custom_expected = intervals(covariance, 7)
  matrix_covariance = Matrix::Matrix(covariance, sparse = TRUE)
  attr(matrix_covariance, "df.t") = 7
  expect_equal(as.numeric(dbreg:::inference_confint(fit, matrix_covariance, c(3, 2), 0.95)),
    as.numeric(custom_expected[c(3, 2), ]), tolerance = 1e-10)
  for (V in list(covariance, unname(covariance))) {
    table = generics::tidy(fit, vcov = V, conf.int = TRUE)
    expect_equal(table[["estimate"]], as.numeric(stats::coef(reference)), tolerance = 1e-10)
    expect_equal(table[["std.error"]], as.numeric(sqrt(diag(covariance))), tolerance = 1e-10)
    expect_equal(table[["conf.low"]], as.numeric(custom_expected[, 1]), tolerance = 1e-10)
    expect_equal(table[["conf.high"]], as.numeric(custom_expected[, 2]), tolerance = 1e-10)
  }
  expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-10)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-10)
  reordered = covariance[c(2, 3, 1), c(2, 3, 1)]
  expect_error(generics::tidy(fit, vcov = reordered), "fitted order")
  wrong = covariance
  dimnames(wrong) = rep(list(c("(Intercept)", "(Intercept)", "grpb")), 2)
  expect_error(generics::tidy(fit, vcov = wrong), "names must match")

  # Unique names still allow a supplied covariance in a different order.
  d[["a"]] = as.numeric(d[["grp"]] == "b")
  unique_fit = dbreg::dbreg(y ~ a + grpb, data = d, strategy = "moments")
  unique_covariance = covariance
  dimnames(unique_covariance) = rep(list(names(stats::coef(unique_fit))), 2)
  unique_covariance = unique_covariance[c(3, 1, 2), c(3, 1, 2)]
  attr(unique_covariance, "df.t") = 7
  table = generics::tidy(unique_fit, vcov = unique_covariance, conf.int = TRUE)
  expect_equal(table[["std.error"]], as.numeric(sqrt(diag(covariance))), tolerance = 1e-10)
  expect_equal(table[["conf.low"]], as.numeric(custom_expected[, 1]), tolerance = 1e-10)

  d[["copy"]] = d[["grpb"]]
  aliased = dbreg::dbreg(y ~ grp + grpb + copy, data = d, strategy = "moments")
  table = generics::tidy(aliased, collin = TRUE, conf.int = TRUE)
  expect_equal(table[["estimate"]][1:3], as.numeric(stats::coef(reference)), tolerance = 1e-10)
  expect_equal(table[["conf.low"]][1:3], as.numeric(expected_iid[, 1]), tolerance = 1e-10)
  expect_true(is.na(table[["estimate"]][4]))
  expect_true(is.na(table[["conf.low"]][4]))

  empty_data = data.frame(y = rep(c(1, 2, 4, 3), 5), g = rep(seq_len(5), each = 4))
  empty_fit = dbreg::dbreg(y ~ 1 | g, data = empty_data, retain = "inference")
  on.exit(close(empty_fit), add = TRUE)
  empty_summary = summary(empty_fit, vcov = "hc1")
  expect_identical(dim(stats::confint(empty_summary)), c(0L, 2L))
  empty_table = generics::tidy(empty_summary, collin = TRUE, conf.int = TRUE)
  expect_equal(nrow(empty_table), 0L)
  expect_identical(empty_table[["conf.low"]], numeric())
  expect_identical(empty_table[["conf.high"]], numeric())

  # An absorbed term between slopes in the formula keeps its own empty row.
  middle = data.frame(x = stats::rnorm(180), z = stats::rnorm(180), g = rep(seq_len(15), each = 12))
  middle[["copy"]] = as.numeric(middle[["g"]])
  middle[["y"]] = 0.5 * middle[["x"]] - 0.3 * middle[["z"]] + cos(middle[["g"]]) + stats::rnorm(180)
  middle_fit = dbreg::dbreg(y ~ x + copy + z | g, data = middle, strategy = "demean")
  middle_reference = stats::lm(y ~ x + z + factor(g), middle)
  table = generics::tidy(middle_fit, collin = TRUE, conf.int = TRUE)
  expect_identical(sort(table[["term"]]), sort(c("x", "copy", "z")))
  positions = match(c("x", "z"), table[["term"]])
  dropped = match("copy", table[["term"]])
  expect_equal(table[["estimate"]][positions], as.numeric(stats::coef(middle_reference)[c("x", "z")]), tolerance = 1e-8)
  expect_equal(table[["conf.low"]][positions], as.numeric(stats::confint(middle_reference)[c("x", "z"), 1]), tolerance = 1e-8)
  expect_true(is.na(table[["estimate"]][dropped]) && is.na(table[["conf.low"]][dropped]))
})
