dbreg_test_sql_prediction_tails = function(conn) {
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources))
  set.seed(929)
  d = data.frame(id = 1:240, x = rep(c(-1, 0, 1), 80), z = sin(1:240), o = 0)
  d[["y"]] = stats::rpois(nrow(d), exp(1 + 2*d[["x"]] - 2*d[["z"]]))
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  fit = dbglm(y ~ x + z, conn = conn, table = table, offset = "o", family = stats::poisson())
  reference = fixest::feglm(y ~ x + z, d, offset = ~o, family = stats::poisson(), notes = FALSE, glm.tol = 1e-12)
  expect_true(all(abs(coef(fit) - coef(reference)) <= 1e-7 + 1e-6*abs(coef(reference))))
  targets = c(-1000, -100, -40, -36, 0, 100, 700, 701, 709, 710)
  data = data.frame(id = seq_along(targets), x = 0, z = 0, o = targets - coef(fit)[["(Intercept)"]])
  # Finite inputs can overflow to Inf-Inf; response means must remain missing.
  overflow = data[1:2, ]
  overflow[["id"]] = length(targets) + 1:2
  overflow[["x"]] = c(1e308, 0)
  overflow[["z"]] = c(1e308, 0)
  overflow[["o"]] = c(0, NA_real_)
  data = rbind(data, overflow)
  inventory = dbglm_validation_tables(conn)
  output = predict(fit, data, conn = conn, row_id = "id", type = "response")
  on.exit(close(output), add = TRUE)
  actual = as.data.frame(output)
  actual = actual[[output[["prediction_column"]]]][match(data[["id"]], actual[["id"]])]
  expected = predict(fit, data, type = "response")
  cancellation = length(targets)+1L
  expect_true(is.na(expected[[cancellation]]))
  expect_true(is.na(predict(fit,data[cancellation,,drop=FALSE],type="response")))
  expect_identical(is.na(actual), is.na(expected))
  expect_identical(is.infinite(actual), is.infinite(expected))
  finite = is.finite(expected)
  expect_true(all(abs(actual[finite]/expected[finite] - 1) < 1e-12))
  expect_identical(actual[1:3], rep(.Machine[["double.eps"]], 3L))
  reference_prediction = as.numeric(predict(reference, data, type = "response"))
  # Inf-Inf is not a reference value: BLAS implementations and batch shape
  # can classify that undefined cancellation differently. Our policy above
  # is missing, with all representable predictor rows still checked to fixest.
  defined = seq_along(actual)!=cancellation
  expect_identical(is.na(actual[defined]), is.na(reference_prediction[defined]))
  expect_identical(is.infinite(actual[defined]), is.infinite(reference_prediction[defined]))
  expect_true(all(abs(actual[finite]/reference_prediction[finite] - 1) < 1e-7))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  missing = DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n FROM ", q(output[["table"]]),
    " WHERE ", q(output[["prediction_column"]]), " IS NULL"))[["n"]]
  expect_equal(as.numeric(missing), 2)
  close(output)
  tail = data[1:2,]
  tail[["x"]] = c(1e308,-1e308)
  tail[["z"]] = 0
  tail[["o"]] = 0
  expected_tail = predict(fit,tail,type="response")
  expect_identical(is.infinite(expected_tail),c(TRUE,FALSE))
  expect_identical(expected_tail[[2L]],.Machine[["double.eps"]])
  tail_output = predict(fit,tail,conn=conn,row_id="id",type="response")
  on.exit(close(tail_output),add=TRUE)
  value = as.data.frame(tail_output)
  value = value[[tail_output[["prediction_column"]]]][match(tail[["id"]],value[["id"]])]
  expect_identical(is.infinite(value),c(TRUE,FALSE))
  expect_identical(value[[2L]],.Machine[["double.eps"]])
  close(tail_output)
  expect_identical(dbglm_validation_tables(conn), inventory)

  # The binomial helper is the existing response kernel, including missing eta.
  eta = c(-1000, -100, -30, -5, 0, 5, 30, 100, 1000, NA_real_)
  frame = data.frame(id = seq_along(eta), prediction = eta)
  input = dbreg:::prediction_upload_frame(conn, frame, resources)
  for (link in c("logit", "cloglog", "probit")) {
    query = dbreg:::sql_prediction_response(conn, gsub('["`]', "", input), "id", "prediction", link)
    actual = dbreg:::backend_get_query(conn, query)
    expected = stats::binomial(link)[["linkinv"]](eta)
    expected = expected[match(actual[["id"]], frame[["id"]])]
    expect_identical(is.na(actual[["prediction"]]), is.na(expected), info = link)
    expect_true(all(abs(actual[["prediction"]] - expected) < 1e-12, na.rm = TRUE), info = link)
  }
}
