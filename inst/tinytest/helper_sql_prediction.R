dbreg_test_sql_prediction = function(conn) {
  before = dbglm_validation_tables(conn)
  original = options(dbreg.prediction_max_entries = 1e6)
  on.exit(options(original), add = TRUE)
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  models = outputs = list()
  on.exit({ for (output in outputs) close(output); for (model in models) close(model) }, add = TRUE)
  set.seed(275)
  n = 480L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
    g = rep(1:12, each = 40), h = sample(1:6, n, TRUE), j = sample(1:3, n, TRUE),
    f = rep(c("a", "O'Brien", "slash\\café"), length.out = n), o = .1*sin(seq_len(n)))
  d[["y"]] = .4*d[["x"]] - .2*d[["z"]] + sin(d[["g"]]) + rnorm(n, sd = .2)
  d[["b"]] = stats::rbinom(n, 1L, stats::plogis(.2*d[["x"]] + .1*sin(d[["g"]])))
  d[["count"]] = stats::rpois(n, exp(.1*d[["x"]] + .1*cos(d[["g"]])))
  table = dbreg:::prediction_upload_frame(conn, d, resources)
  check = function(fit, data, type = "link", sql_source = FALSE) {
    args = list(object = fit, conn = conn, row_id = "id")
    if (inherits(fit, "dbglm")) args[["type"]] = type
    if (sql_source) args[["table"]] = dbreg:::prediction_upload_frame(conn, data, resources) else args[["newdata"]] = data
    output = do.call(predict, args)
    actual = as.data.frame(output)
    reference = if (inherits(fit, "dbglm")) predict(fit, data, type = type) else predict(fit, data)
    reference = reference[match(actual[["id"]], data[["id"]])]
    expect_identical(is.na(actual[[output[["prediction_column"]]]]), is.na(reference))
    expect_true(all(abs(actual[[output[["prediction_column"]]]] - reference) <= 1e-8 + 1e-6*abs(reference), na.rm = TRUE))
    expect_equal(nrow(actual), nrow(data))
    expect_equal(nrow(head(output, 3)), min(3L, nrow(data)))
    expect_equal(nrow(head(output, 0)), 0L)
    outputs[[length(outputs) + 1L]] = output
    close(output)
    expect_error(head(output), "released")
    invisible(NULL)
  }
  message("[SQL prediction] fitted factors, powers and interactions")
  fit = dbreg(y ~ x + f + I(x^2) + x:f, conn = conn, table = table, strategy = "moments")
  models[[length(models) + 1L]] = fit
  reference_data = d
  reference_data[["f"]] = factor(d[["f"]], levels = fit[["factor_levels"]][["f"]])
  reference = fixest::feols(y ~ x + f + I(x^2) + x:f, reference_data, notes = FALSE)
  # fixest 0.14.0 adds a redundant I() wrapper around power terms.
  expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-8)
  newdata = d[41:57, setdiff(names(d), c("y", "b", "count")), drop = FALSE]
  newdata[["x"]][2] = NA_real_
  check(fit, newdata)
  check(fit, newdata, sql_source = TRUE)
  decimal = d
  decimal[["category"]] = rep(c(.12345678901234567, .76543210987654321, 1.2345678901234567), length.out = n)
  decimal_table = dbreg:::prediction_upload_frame(conn, decimal, resources)
  numeric_factor = dbreg(y ~ x + category, conn = conn, table = decimal_table,
    factor_refs = list(category = decimal[["category"]][[1L]]), strategy = "moments")
  models[[length(models) + 1L]] = numeric_factor
  check(numeric_factor, decimal[1:12, ])
  check(numeric_factor, decimal[1:12, ], sql_source = TRUE)
  decimal[["category"]] = factor(decimal[["category"]], levels = numeric_factor[["factor_levels"]][["category"]])
  numeric_reference = stats::lm(y ~ x + category, decimal)
  expect_equal(unname(coef(numeric_factor)), unname(coef(numeric_reference)), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(numeric_factor)), as.numeric(vcov(numeric_reference)), tolerance = 1e-8)
  bad = newdata
  bad[["f"]][3] = "new category"
  inventory = dbglm_validation_tables(conn)
  expect_error(predict(fit, bad, conn = conn, row_id = "id"), "unseen factor")
  expect_identical(dbglm_validation_tables(conn), inventory)
  bad = newdata
  bad[["id"]][3] = bad[["id"]][1]
  expect_error(predict(fit, bad, conn = conn, row_id = "id"), "unique and nonmissing")
  expect_identical(dbglm_validation_tables(conn), inventory)
  bad[["id"]][3] = NA_integer_
  expect_error(predict(fit, bad, conn = conn, row_id = "id"), "unique and nonmissing")
  expect_identical(dbglm_validation_tables(conn), inventory)
  expect_error(predict(fit, newdata, conn = conn), "explicit.*row_id")
  expect_error(predict(fit, newdata, conn = conn, row_id = "id", interval = "confidence"), "intervals")
  empty = newdata[FALSE, , drop = FALSE]
  check(fit, empty)
  export = basename(tempfile("prediction_export_"))
  out = predict(fit, newdata, conn = conn, row_id = "id", prediction_table = export)
  close(out)
  expect_true(dbreg:::sample_table_exists(conn, export))
  expect_error(predict(fit, newdata, conn = conn, row_id = "id", prediction_table = export), "already exists")
  dbglm_validation_drop(conn, export)

  message("[SQL prediction] exact BIGINT/date tuples and nonstandard output keys")
  typed = basename(tempfile("prediction_typed_"))
  resources[["tables"]] = c(resources[["tables"]], typed)
  quote_id = function(value) as.character(DBI::dbQuoteIdentifier(conn, value))
  typed_sql = paste0("SELECT *,CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS keyg,",
    "CAST('2000-01-01' AS DATE)+CAST(g AS INTEGER) AS datekey,",
    "CAST(9007199254740992 AS BIGINT)+CAST(id AS BIGINT) AS ", quote_id("row key"), " FROM ", table)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, typed, typed_sql, c(names(d), "keyg", "datekey", "row key")))
  fit = dbreg(y ~ x | keyg^datekey, conn = conn, table = quote_id(typed), strategy = "demean", retain = "prediction")
  models[[length(models) + 1L]] = fit
  typed_data = d[1:8, c("id", "x", "g")]
  typed_data[["keyg"]] = bit64::as.integer64("9007199254740992") + bit64::as.integer64(typed_data[["g"]])
  typed_data[["datekey"]] = as.Date("2000-01-01") + typed_data[["g"]]
  typed_data[["row key"]] = bit64::as.integer64("9007199254740992") + bit64::as.integer64(typed_data[["id"]])
  output = predict(fit, typed_data, conn = conn, row_id = "row key")
  actual = as.data.frame(output)
  index = match(as.character(actual[["row key"]]), as.character(typed_data[["row key"]]))
  expect_false(anyNA(index))
  expect_equal(actual[[output[["prediction_column"]]]], predict(fit, typed_data)[index], tolerance = 1e-8)
  key_type = DBI::dbGetQuery(conn, paste("SELECT typeof(", quote_id("row key"), ") AS type FROM", quote_id(output[["table"]]), "LIMIT 1"))[[1]][[1L]]
  expect_identical(tolower(key_type), "bigint")
  close(output)
  typed_data[["keyg"]] = as.character(typed_data[["keyg"]])
  typed_data[["datekey"]] = as.character(typed_data[["datekey"]])
  typed_data[["keyg"]][2] = "9.007199254740993e15"
  typed_data[["datekey"]][3] = paste0(typed_data[["datekey"]][3], " 12:00:00")
  output = predict(fit, typed_data, conn = conn, row_id = "id")
  actual = as.data.frame(output)
  expect_identical(is.na(actual[[output[["prediction_column"]]]]), is.na(predict(fit, typed_data)[match(actual[["id"]], typed_data[["id"]])]))
  expect_equal(sum(is.na(actual[[output[["prediction_column"]]]])), 2L)
  close(output)
  output = predict(fit, conn = conn, table = quote_id(typed), row_id = "row key")
  actual = as.data.frame(output)
  expect_equal(nrow(actual), n)
  expect_identical(sort(as.character(actual[["row key"]])), sort(as.character(bit64::as.integer64("9007199254740992") + bit64::as.integer64(d[["id"]]))))
  close(output)

  message("[SQL prediction] portable multi-FE effects and varying slopes")
  for (fml in list(y ~ x + z | g + h + j)) {
    fit = dbreg(fml, conn = conn, table = table, strategy = "demean", retain = "prediction")
    models[[length(models) + 1L]] = fit
    changed = newdata
    changed[["g"]][3] = 999L
    check(fit, changed)
    check(fit, changed, sql_source = TRUE)
  }
  message("[SQL prediction] SQL dictionaries and owned output lifetimes")
  options(dbreg.prediction_max_entries = 64)
  fit = dbreg(y ~ x + z | g + h, conn = conn, table = table, strategy = "demean", retain = "prediction")
  models[[length(models) + 1L]] = fit
  expect_identical(fit[["fitted_state"]][["prediction"]][["storage"]], "database")
  check(fit, newdata)
  out = predict(fit, newdata, row_id = "id")
  outputs[[length(outputs) + 1L]] = out
  expected = as.data.frame(out)
  close(fit)
  expect_equal(as.data.frame(out), expected)
  restored = unserialize(serialize(out, NULL))
  expect_error(head(restored), "live owning session")
  close(restored)
  expect_equal(as.data.frame(out), expected)
  close(out)
  expect_true(DBI::dbIsValid(conn))
  options(dbreg.prediction_max_entries = 1e6)

  for (family in list(stats::binomial("logit"), stats::binomial("cloglog"), stats::binomial("probit"), stats::poisson())) {
    message("[SQL prediction] ", family[["family"]], "/", family[["link"]])
    response = if (family[["family"]] == "poisson") "count" else "b"
    fit = dbglm(stats::as.formula(paste0(response, "~x+f|g")), conn = conn, table = table,
      family = family, offset = "o", fixef_rm = "none", retain = "prediction",
      control = dbglm_control(score_tol = 1e-8, glm_tol = 1e-10))
    models[[length(models) + 1L]] = fit
    reference_data = d
    reference_data[["f"]] = factor(d[["f"]], levels = fit[["model_spec"]][["factor_levels"]][["f"]])
    reference = fixest::feglm(stats::as.formula(paste0(response, "~x+f|g")), reference_data,
      family = family, offset = ~o, fixef.rm = "none", notes = FALSE, glm.tol = 1e-12, fixef.tol = 1e-10)
    expect_true(all(abs(unname(coef(fit)) - unname(coef(reference))) <= 1e-7 + 1e-6*abs(coef(reference))),
      info = paste(family[["link"]], "max coefficient error", max(abs(unname(coef(fit)) - unname(coef(reference))))))
    expect_true(all(abs(as.numeric(vcov(fit)) - as.numeric(vcov(reference))) <= 1e-8 + 1e-5*abs(as.numeric(vcov(reference)))), info = family[["link"]])
    changed = newdata
    changed[["g"]][3] = 999L
    changed[["o"]][4:5] = c(-100, 100)
    check(fit, changed, "link")
    check(fit, changed, "response", sql_source = TRUE)
  }
  fit = dbglm(count ~ x | g[z] + h, conn = conn, table = table,
    family = stats::poisson(), fixef_rm = "none", retain = "prediction")
  models[[length(models) + 1L]] = fit
  check(fit, newdata)
  for (output in outputs) close(output)
  for (model in models) close(model)
  dbreg:::release_fit_resources(resources)
  resources[["tables"]] = character()
  expect_identical(dbglm_validation_tables(conn), before)
}
