# In-memory BIGINT keys must remain exact, including adjacent values above 2^53.
# DuckDB's numeric conversion mode can otherwise import the stored bits as tiny
# doubles. Default connections preserve the type; supplied connections either
# preserve it or fail before estimation, without leaking registered relations.
if (requireNamespace("bit64", quietly = TRUE)) local({
  set.seed(7352)
  n = 600L
  groups = rep(1:6, length.out = n)
  labels = bit64::as.integer64(c("9007199254740992", "9007199254740993",
    "9007199254740994", "9007199254740995", "9007199254740996", "9007199254740997"))
  d = data.frame(id = bit64::as.integer64("9007199254740992") + bit64::as.integer64(seq_len(n) - 1L), g = labels[groups],
    x = rnorm(n), z = rnorm(n), a = bit64::as.integer64(rep(1:3, length.out = n)))
  d[["y"]] = .2 + .5 * d[["x"]] - .3 * d[["z"]] + sin(groups) + rnorm(n)
  d[["binary"]] = rbinom(n, 1, binomial("cloglog")[["linkinv"]](
    -1 + .3 * d[["x"]] - .2 * d[["z"]] + .2 * sin(groups)))
  reference_data = d
  reference_data[["g"]] = as.character(d[["g"]])
  reference_data[["a"]] = as.numeric(d[["a"]])
  reference = fixest::feols(y ~ x + z | g, reference_data, weights = ~a, fixef.rm = "none")
  reference_glm = fixest::feglm(binary ~ x + z | g, reference_data,
    weights = ~a, family = binomial("cloglog"), fixef.rm = "none", glm.tol = 1e-12)
  check = function(conn = NULL) {
    fit = dbreg(y ~ x + z | g, data = d, conn = conn, weights = "a", strategy = "demean")
    expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-8)
    export = if (!is.null(conn)) "integer64_sample" else NULL
    fit_glm = dbglm(binary ~ x + z | g, data = d, conn = conn, weights = "a", row_id = "id", sample_table = export)
    expect_equal(coef(fit_glm), coef(reference_glm), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit_glm)), as.numeric(vcov(reference_glm)), tolerance = 1e-7)
    expect_equal(fit_glm[["n_params_model"]], reference_glm[["nparams"]])
    if (!is.null(export)) {
      keys = DBI::dbGetQuery(conn, paste("SELECT CAST(id AS VARCHAR) AS id FROM", export, "ORDER BY id"))[["id"]]
      expect_identical(keys, as.character(d[["id"]]))
      DBI::dbRemoveTable(conn, export)
    }
  }
  check()
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "integer64", config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  check(conn)
  expect_equal(length(DBI::dbListTables(conn)), 0L)
  # A supplied connection may support integer64 independently of its output
  # mode in future drivers. Compare results in that case; otherwise fail closed.
  numeric_conn = DBI::dbConnect(duckdb::duckdb(bigint = "numeric", config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(numeric_conn, shutdown = TRUE), add = TRUE)
  duckdb::duckdb_register(numeric_conn, "type_probe", d)
  schema = DBI::dbGetQuery(numeric_conn, "DESCRIBE type_probe")
  preserves = identical(schema[["column_type"]][match("g", schema[["column_name"]])], "BIGINT")
  duckdb::duckdb_unregister(numeric_conn, "type_probe")
  if (preserves) check(numeric_conn) else {
    expect_error(dbreg(y ~ x | g, data = d, conn = numeric_conn), "did not preserve integer64")
    expect_error(dbglm(binary ~ x | g, data = d, conn = numeric_conn), "did not preserve integer64")
  }
  expect_true(DBI::dbIsValid(numeric_conn))
  expect_equal(length(DBI::dbListTables(numeric_conn)), 0L)
  # Prediction must decode numeric integer64 values while preserving exact
  # factor labels. Negative integer64 storage otherwise appears as NaN to R.
  d[["q"]] = bit64::as.integer64(rep(-3:3, length.out = n))
  reference_data[["q"]] = as.numeric(d[["q"]])
  d[["g_factor"]] = factor(as.character(d[["g"]]))
  reference_data[["g_factor"]] = d[["g_factor"]]
  linear = dbreg(y ~ q + g_factor, data = d, strategy = "moments")
  linear_ref = lm(y ~ q + g_factor, reference_data)
  expect_equal(unname(predict(linear, d)), unname(predict(linear_ref, reference_data)), tolerance = 1e-8)
  binary = dbglm(binary ~ q + factor(g), data = d)
  binary_ref = glm(binary ~ q + factor(g), reference_data, family = binomial("cloglog"),
    control = glm.control(epsilon = 1e-12))
  expect_equal(unname(predict(binary, d, type = "link")), unname(predict(binary_ref, reference_data, type = "link")), tolerance = 1e-7)
  expect_equal(unname(predict(binary, d, type = "response")), unname(predict(binary_ref, reference_data, type = "response")), tolerance = 1e-7)
  missing = d
  missing[["q"]][1] = NA
  expected_missing = predict(binary_ref, reference_data, type = "link")
  expected_missing[1] = NA
  expect_equal(unname(predict(binary, missing)), unname(expected_missing), tolerance = 1e-7)
  d[["o"]] = bit64::as.integer64(rep(-1:1, length.out = n))
  reference_data[["o"]] = as.numeric(d[["o"]])
  offset_fit = dbglm(binary ~ q, data = d, offset = "o")
  offset_ref = glm(binary ~ q, reference_data, offset = reference_data[["o"]],
    family = binomial("cloglog"), control = glm.control(epsilon = 1e-12))
  expect_equal(unname(predict(offset_fit, d)), unname(predict(offset_ref, reference_data, type = "link")), tolerance = 1e-7)
  d[["count"]] = bit64::as.integer64(round(20 + d[["y"]]))
  reference_data[["count"]] = as.numeric(d[["count"]])
  within_fit = dbreg(count ~ x | g, data = d, strategy = "demean")
  within_ref = lm(count ~ x + g_factor, reference_data)
  expect_equal(unname(predict(within_fit, d)), unname(predict(within_ref, reference_data)), tolerance = 1e-8)
})
