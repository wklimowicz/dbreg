dbreg_test_sql_prediction_fe_only = function(conn) {
  set.seed(1827)
  n = 240L
  d = data.frame(id = seq_len(n), g = rep(1:8, length.out = n), h = sample.int(4L, n, TRUE), j = sample.int(3L, n, TRUE))
  d[["y"]] = sin(d[["g"]])+cos(d[["h"]])+rnorm(n)
  table = dbglm_validation_upload(conn, d)
  on.exit(dbglm_validation_drop(conn, table))
  for (formula in list(y ~ 1 | g + h, y ~ 1 | g + h + j)) {
    fit = dbreg::dbreg(formula, retain = "all", conn = conn, table = table, fixef_rm = "none")
    reference = fixest::feols(formula, d, fixef.rm = "none")
    expect_equal(predict(fit, d), predict(reference, d), tolerance = 1e-7)
    output = predict(fit, newdata = d[1:17, ], conn = conn, row_id = "id")
    actual = as.data.frame(output)
    expected = predict(reference, d)[match(actual[["id"]], d[["id"]])]
    expect_equal(actual[[output[["prediction_column"]]]], expected, tolerance = 1e-7)
    expect_identical(anyNA(actual[[output[["prediction_column"]]]]), FALSE)
    close(output)
    close(fit)
  }
}
