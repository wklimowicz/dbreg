# Value-based checks share the same synthetic rows on DuckDB and Databricks.
dbreg_test_linear_prediction = function(conn) {
  before = dbglm_validation_tables(conn)
  set.seed(819)
  n = 480L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
    g = rep(1:8, each = 60), h = rep(1:5, length.out = n),
    k = rep(1:3, each = 160), f = rep(c("a", "b", "c"), 160),
    w = rep(c(.25, .5, 1), 160))
  d[["y"]] = 3 + .6 * d[["x"]] - .2 * d[["z"]] + .3 * (d[["f"]] == "b") +
    sin(d[["g"]]) + .2 * cos(d[["h"]]) + rnorm(n, sd = .2)
  uploaded = dbglm_validation_upload(conn, d)
  table = basename(tempfile("linear_prediction_"))
  on.exit({dbglm_validation_drop(conn, table); dbglm_validation_drop(conn, uploaded)})
  fields = names(d)
  fields[fields == "x"] = "CASE WHEN id=2 THEN CAST(NULL AS DOUBLE) ELSE x END AS x"
  fields[fields == "y"] = "CASE WHEN id=3 THEN CAST(NULL AS DOUBLE) ELSE y END AS y"
  fields[fields == "w"] = "CASE WHEN id=4 THEN 0.0 ELSE w END AS w"
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table,
    paste("SELECT", paste(fields, collapse = ","), "FROM", uploaded), names(d)))
  d[["x"]][2] = NA_real_
  d[["y"]][3] = NA_real_
  d[["w"]][4] = 0
  reference_data = d
  reference_data[["f"]] = factor(d[["f"]], levels = c("c", "a", "b"))
  inventory = sort(c(before, uploaded, table))
  check = function(fml, strategy) {
    message("[OLS prediction] ", strategy, ": ", paste(deparse(fml), collapse = ""))
    fit = dbreg(fml, conn = conn, table = table, weights = "w", strategy = strategy,
      factor_refs = list(f = "c"), retain = "prediction")
    compact = dbreg(fml, conn = conn, table = table, weights = "w", strategy = strategy, factor_refs = list(f = "c"))
    ref = fixest::feols(fml, reference_data, weights = ~w, fixef.rm = "none", notes = FALSE)
    expect_equal(coef(fit), coef(compact), tolerance = 1e-7)
    expect_equal(vcov(fit), vcov(compact), tolerance = 1e-7)
    reference_coef = coef(ref)
    # fixest 0.14.0 redundantly wraps I() term labels; keep native dbreg names.
    names(reference_coef) = sub("^I\\(I\\((.*)\\)\\)$", "I(\\1)", names(reference_coef))
    expect_equal(coef(fit), reference_coef, tolerance = 1e-6)
    expect_equal(nobs(fit), nobs(ref))
    newdata = d[rev(seq_len(n)), setdiff(names(d), c("y", "w")), drop = FALSE]
    found = predict(fit, newdata)
    expected = unname(predict(ref, newdata))
    expect_equal(found, expected, tolerance = 1e-6, info = paste(strategy, deparse(fml)))
    expect_equal(predict(fit, newdata[1, , drop = FALSE]), found[1], tolerance = 1e-10)
    expect_equal(predict(fit, newdata[c(3, 1, 3), , drop = FALSE]), found[c(3, 1, 3)], tolerance = 1e-10)
    expect_equal(predict(fit, newdata[FALSE, , drop = FALSE]), numeric())
    changed = newdata
    changed[["y"]] = 1e9
    changed[["w"]] = 0
    expect_equal(predict(fit, changed), found, tolerance = 1e-12)
    expect_equal(predict(unserialize(serialize(fit, NULL)), newdata), found, tolerance = 1e-12)
    # A missing fitting outcome and a zero fitting weight still allow prediction.
    positions = match(c(2, 3, 4), newdata[["id"]])
    expect_true(is.na(found[positions[1]]) && all(is.finite(found[positions[2:3]])))
    if (fit[["strategy"]] == "compress") {
      legacy = fit
      legacy[["fitted_state"]] = NULL
      expect_equal(predict(fit), predict(legacy), tolerance = 1e-12)
      expect_equal(predict(fit, fit[["data"]]), predict(legacy), tolerance = 1e-12)
    } else expect_error(predict(fit), "newdata")
    expect_error(predict(fit, newdata, interval = "confidence"), "Intervals")
    expect_identical(dbglm_validation_tables(conn), inventory)
    fit
  }
  for (strategy in c("demean", "compress", "auto")) {
    fit = check(y ~ x + z + f + I(x^2) | g + h, strategy)
  }
  check(y ~ x:f + z | g + h, "demean")
  check(y ~ x + f | g^h, "compress")
  check(y ~ x + f | g + h + k, "demean")
  # Prediction must not use a later version of the caller's source.
  input = d[1:8, setdiff(names(d), c("y", "w")), drop = FALSE]
  saved = predict(fit, input)
  dbglm_validation_drop(conn, table)
  expect_equal(predict(fit, input), saved, tolerance = 1e-12)
  expect_identical(dbglm_validation_tables(conn), sort(c(before, uploaded)))
  expect_true(DBI::dbIsValid(conn))
}

dbreg_test_linear_prediction_keys = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(g = rep(1:4, each = 40), x = rep(c(-1, 0, 1, 2), 40))
  d[["y"]] = 2 + .6 * d[["x"]] + d[["g"]] + sin(seq_len(nrow(d))) / 10
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("linear_prediction_keys_"))
  on.exit({dbglm_validation_drop(conn, table); dbglm_validation_drop(conn, source)})
  sql = paste0("SELECT *,CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS big,",
    "CAST(1.0 AS DOUBLE)+g*CAST(2.220446049250313e-16 AS DOUBLE) AS close_key,",
    "CASE WHEN g IN (1,2) THEN 'a|b' ELSE 'a' END AS part FROM ", source)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, sql, c(names(d), "big", "close_key", "part")))
  newdata = d[c(1, 41, 81, 121), c("g", "x")]
  newdata[["big"]] = bit64::as.integer64("9007199254740992") + bit64::as.integer64(newdata[["g"]])
  newdata[["close_key"]] = 1 + newdata[["g"]] * 2.220446049250313e-16
  newdata[["part"]] = ifelse(newdata[["g"]] %in% c(1, 2), "a|b", "a")
  ref = fixest::feols(y ~ x | g, d, fixef.rm = "none", notes = FALSE)
  expected = unname(predict(ref, newdata))
  for (strategy in c("demean", "compress")) for (fml in list(y ~ x | big, y ~ x | close_key^part)) {
    message("[OLS prediction keys] ", strategy, ": ", paste(deparse(fml), collapse = ""))
    fit = dbreg(fml, conn = conn, table = table, strategy = strategy, retain = "prediction")
    compact = dbreg(fml, conn = conn, table = table, strategy = strategy)
    expect_equal(coef(compact), coef(ref), tolerance = 1e-7)
    expect_true(all(abs(as.numeric(vcov(compact)) - as.numeric(vcov(ref))) <= 1e-10 + 1e-5 * abs(as.numeric(vcov(ref)))))
    expect_identical(dimnames(vcov(compact)), dimnames(vcov(ref)))
    expect_equal(compact[["df_residual"]], df.residual(ref))
    expect_equal(attr(compact[["vcov"]], "rss"), sum(residuals(ref)^2), tolerance = 1e-7)
    expect_equal(predict(fit, newdata), expected, tolerance = 1e-6)
    expect_equal(predict(fit, newdata[1, , drop = FALSE]), expected[1], tolerance = 1e-6)
    expect_identical(length(fit[["fitted_state"]][["specification"]][["dictionaries"]][[1]][["key"]]), 4L)
  }
  expect_identical(dbglm_validation_tables(conn), sort(c(before, source, table)))
}

dbreg_test_linear_prediction_estimability = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(g = rep(1:4, each = 40), h = rep(1:4, each = 40),
    x = rep(c(-1, 0, 1, 2), 40))
  d[["z"]] = 2 * d[["x"]]
  d[["f"]] = c("a", "b", "c", "d")[d[["g"]]]
  d[["y"]] = 2 + .6 * d[["x"]] + d[["g"]] + sin(seq_len(nrow(d))) / 10
  table = dbglm_validation_upload(conn, d)
  on.exit(dbglm_validation_drop(conn, table))
  for (strategy in c("demean", "compress")) {
    message("[OLS estimability] ", strategy)
    fit = dbreg(y ~ x + z | g + h, conn = conn, table = table, strategy = strategy, retain = "prediction")
    rows = d[c(1, 41), c("x", "z", "g", "h")]
    expect_true(all(is.finite(predict(fit, rows))))
    crossed = rows
    crossed[["h"]] = rev(crossed[["h"]])
    expect_true(all(is.na(predict(fit, crossed))))
    rows[["z"]][2] = rows[["z"]][2] + 1
    expect_true(is.finite(predict(fit, rows)[1]) && is.na(predict(fit, rows)[2]))
    rows[["g"]][1] = 99
    expect_true(all(is.na(predict(fit, rows))))
  }
  # An FE key can also appear as a structurally absorbed ordinary factor term.
  fit = dbreg(y ~ x + f | f, conn = conn, table = table, strategy = "demean", retain = "prediction")
  rows = d[1:2, c("x", "f")]
  rows[["f"]][2] = "new"
  expect_true(is.finite(predict(fit, rows)[1]) && is.na(predict(fit, rows)[2]))
  expect_error(dbreg(y ~ x | g, conn = conn, table = table, retain = "prediction", strategy = "mundlak"), "Mundlak")
  expect_error(dbreg(y ~ x | g, conn = conn, table = table, retain = "prediction", sql_only = TRUE), "sql_only")
  expect_error(suppressWarnings(dbreg(y ~ x | g, conn = conn, table = table, retain = "prediction", query_only = TRUE)), "sql_only")
  expect_error(dbreg(y ~ x | g, conn = conn, table = table, retain = "prediction", data_only = TRUE), "data_only")
  original = getOption("dbreg.prediction_max_entries")
  on.exit(options(dbreg.prediction_max_entries = original), add = TRUE)
  options(dbreg.prediction_max_entries = 1)
  expect_error(dbreg(y ~ x | g, conn = conn, table = table, retain = "prediction", strategy = "demean"), "prediction state exceeds")
  options(dbreg.prediction_max_entries = original)
  expect_identical(dbglm_validation_tables(conn), sort(c(before, table)))
}

dbreg_test_cluster_key_precision = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(g = rep(1:4, each = 40), x = rep(c(-1, 0, 1, 2), 40), cl = rep(1:8, 20))
  d[["y"]] = 2 + .6 * d[["x"]] + d[["g"]] + sin(seq_len(nrow(d))) / 10
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("linear_cluster_keys_"))
  on.exit({dbglm_validation_drop(conn, table); dbglm_validation_drop(conn, source)})
  sql = paste0("SELECT *,CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS big,",
    "CAST(9007199254740992 AS BIGINT)+CAST(cl AS BIGINT) AS big_cluster,",
    "CAST(1.0 AS DOUBLE)+cl*CAST(2.220446049250313e-16 AS DOUBLE) AS close_cluster FROM ", source)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, sql, c(names(d), "big", "big_cluster", "close_cluster")))
  reference = fixest::feols(y ~ x | g, d, vcov = ~cl, fixef.rm = "none", notes = FALSE,
    ssc = fixest::ssc(K.fixef = "full"))
  for (cluster in c("big_cluster", "close_cluster")) {
    message("[OLS cluster keys] ", cluster)
    fit = dbreg(y ~ x | big, conn = conn, table = table, strategy = "compress", vcov = stats::reformulate(cluster))
    expect_equal(coef(fit), coef(reference), tolerance = 1e-7)
    expect_true(all(abs(as.numeric(vcov(fit)) - as.numeric(vcov(reference))) <= 1e-10 + 1e-5 * abs(as.numeric(vcov(reference)))))
    expect_equal(fit[["n_clusters"]], 8)
  }
  robust = dbreg(y ~ x | big, conn = conn, table = table, strategy = "compress", vcov = "hc1")
  reference = fixest::feols(y ~ x | g, d, vcov = "HC1", fixef.rm = "none", notes = FALSE,
    ssc = fixest::ssc(K.fixef = "full"))
  expect_true(all(abs(as.numeric(vcov(robust)) - as.numeric(vcov(reference))) <= 1e-10 + 1e-5 * abs(as.numeric(vcov(reference)))))
  expect_identical(dbglm_validation_tables(conn), sort(c(before, source, table)))
}

dbreg_test_linear_prediction_snapshot = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(g = rep(1:4, each = 40), x = rep(c(-1, 0, 1, 2), 40), f = rep(c("a", "b"), 80))
  d[["y"]] = 2 + .6 * d[["x"]] + .2 * (d[["f"]] == "b") + d[["g"]] + sin(seq_len(nrow(d))) / 10
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("linear_prediction_snapshot_"))
  on.exit({dbglm_validation_drop(conn, table); dbglm_validation_drop(conn, source)}, add = TRUE)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  columns = c("group label", "x value", "f", "y")
  sql = paste("SELECT", paste(paste(c("g", "x", "f", "y"), "AS", q(columns)), collapse = ","), "FROM", source)
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, sql, columns))
  fml = y ~ `x value` + f | `group label`
  compact = dbreg(fml, conn = conn, table = table, strategy = "demean")
  original = get("finalize_linear_prediction", asNamespace("dbreg"))
  recover = function(result, inputs) {
    dbglm_validation_drop(conn, table)
    original(result, inputs)
  }
  assignInNamespace("finalize_linear_prediction", recover, "dbreg")
  on.exit(assignInNamespace("finalize_linear_prediction", original, "dbreg"), add = TRUE)
  fit = dbreg(fml, conn = conn, table = table, strategy = "demean", retain = "prediction")
  assignInNamespace("finalize_linear_prediction", original, "dbreg")
  names(d) = columns
  newdata = d[1:8, setdiff(columns, "y"), drop = FALSE]
  ref = fixest::feols(fml, d, fixef.rm = "none", notes = FALSE)
  expect_equal(coef(fit), coef(compact), tolerance = 1e-7)
  expect_equal(predict(fit, newdata), unname(predict(ref, newdata)), tolerance = 1e-6)
  expect_identical(dbglm_validation_tables(conn), sort(c(before, source)))
  expect_true(DBI::dbIsValid(conn))
}
