# GLM source binding is independent of the physical SQL/design namespace.
# Explicit raw X and glm.fit supply the statistical oracle.
local({
  set.seed(29013)
  proposed = sprintf("dbreg_internal_x_%04d__", seq_len(3L))
  expect_identical(dbreg:::allocate_sql_aliases(proposed, rep(TRUE, 3L)), proposed)
  reserved = toupper(c(proposed[1:2], paste0(proposed[[1L]], "_1")))
  allocated = dbreg:::allocate_sql_aliases(proposed, rep(TRUE, 3L), reserved_names = reserved)
  expect_false(any(tolower(allocated) %in% tolower(reserved)))
  expect_identical(anyDuplicated(tolower(allocated)), 0L)
  expect_identical(allocated[[3L]], proposed[[3L]])
  n = 320L
  d = data.frame(id = seq_len(n), check.names = FALSE)
  d[["x value"]] = rnorm(n, sd = .4)
  d[["x_value"]] = rnorm(n, sd = .4)
  d[["group name"]] = sample(c("a_x_b", "a:b", "c"), n, TRUE)
  d[["numeric group"]] = sample(1:3, n, TRUE)
  d[["period value"]] = rep(1:4, length.out = n)
  d[["treat value"]] = rbinom(n, 1, .5)
  d[["weight value"]] = rep(c(.5, 1, 2, 3), length.out = n)
  d[["off value"]] = runif(n, -.2, .2)
  d[["cluster value"]] = rep(seq_len(16L), each = 20L)
  d[["future cluster"]] = rep(seq_len(10L), length.out = n)
  d[["in sample"]] = seq_len(n) %% 7L != 0L
  levels = c("a:b", "a_x_b", "c")
  X = cbind(1, d[["x value"]], d[["x_value"]], d[["x value"]]^2,
    sapply(levels[-1L], function(level) as.numeric(d[["group name"]] == level)),
    as.numeric(d[["numeric group"]] == 2), as.numeric(d[["numeric group"]] == 3),
    sapply(2:4, function(period) as.numeric(d[["period value"]] == period) * d[["treat value"]]))
  eta = as.numeric(X %*% seq(.1, .3, length.out = ncol(X))) + d[["off value"]]
  fml = `y value` ~ `x value` + x_value + I(`x value`^2) + `group name` +
    factor(`numeric group`) + i(`period value`, `treat value`, ref = 1)
  bare_ssc = list(K.adj = FALSE, G.adj = FALSE)
  covariance = function(X, y, w, mu, family, cluster = NULL, iid = FALSE) {
    derivative = family[["mu.eta"]](family[["linkfun"]](mu))
    variance = family[["variance"]](mu)
    bread = solve(crossprod(X, (w * derivative^2 / variance) * X))
    if (iid) return(bread * if (family[["family"]] == "gaussian")
      sum(w * (y - mu)^2) / (length(y) - ncol(X)) else 1)
    scores = X * (w * derivative / variance * (y - mu))
    if (!is.null(cluster)) scores = rowsum(scores, cluster, reorder = FALSE)
    bread %*% crossprod(scores) %*% bread
  }
  minimal_names = c("id", "x value", "x_value", "group name", "numeric group",
    "period value", "treat value", "off value")
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  for (family in list(gaussian(), binomial(), poisson())) {
    d[["y value"]] = if (family[["family"]] == "gaussian") eta + rnorm(n, sd = .5) else
      if (family[["family"]] == "binomial") rbinom(n, 1, plogis(eta)) else rpois(n, exp(eta))
    eligible = d[["in sample"]]
    reference = suppressWarnings(glm.fit(X[eligible,,drop = FALSE], d[["y value"]][eligible],
      weights = d[["weight value"]][eligible], offset = d[["off value"]][eligible],
      family = family, control = glm.control(epsilon = 1e-12, maxit = 100)))
    expected_eta = as.numeric(X %*% reference[["coefficients"]]) + d[["off value"]]
    for (strategy in c("moments", "compress")) for (kind in c("iid", "hc1", "cluster")) {
      fit = dbglm(fml, data = d, family = family, weights = ~`weight value`, offset = ~`off value`,
        subset = ~`in sample`, vcov = if (kind == "cluster") ~`cluster value` else kind,
        ssc = bare_ssc, strategy = strategy, control = dbglm_control(glm_tol = 1e-10, score_tol = 1e-9))
      expect_equal(unname(coef(fit)), unname(reference[["coefficients"]]), tolerance = 1e-6)
      expected_V = covariance(X[eligible,,drop = FALSE], d[["y value"]][eligible],
        d[["weight value"]][eligible], reference[["fitted.values"]], family,
        cluster = if (kind == "cluster") d[["cluster value"]][eligible], iid = kind == "iid")
      expect_equal(as.numeric(vcov(fit)), as.numeric(expected_V), tolerance = 1e-6)
      expect_equal(predict(fit, d, type = "link"), expected_eta, tolerance = 1e-6)
      minimal = d[minimal_names]
      original = minimal
      expect_equal(predict(fit, minimal, type = "link"), expected_eta, tolerance = 1e-6)
      expect_identical(minimal, original)
      expect_identical(fit[["model_spec"]][["source_vars"]], setdiff(minimal_names, c("id", "off value")))
      expect_false(any(c("y value", "weight value", "cluster value", "in sample") %in%
        dbreg:::prediction_required_columns(fit)))
      expect_identical(attr(vcov(fit), "dbreg_coefficient_ids"), fit[["model_spec"]][["fitted_ids"]])
      expect_identical(fit[["linear_design"]], NULL)
      expect_identical(fit[["model_spec"]][["columns"]][["alias"]],
        sprintf("dbreg_internal_x_%04d__", seq_len(nrow(fit[["model_spec"]][["columns"]]))))
      aliases = fit[["model_spec"]][["source_map"]]
      # x_value is a separate raw source, so a generated alias must not alter
      # either original input. Add extras from an untouched raw frame.
      extra = original
      physical = names(aliases)[aliases == "x value"]
      for (name in setdiff(physical, names(original))) extra[[name]] = 999
      expect_equal(predict(fit, extra, type = "link"), expected_eta, tolerance = 1e-6)
      absent = extra
      absent[["x value"]] = NULL
      expect_error(predict(fit, absent), "x value")
      unseen = original[1:3,]
      unseen[["group name"]] = "unseen"
      expect_error(predict(fit, unseen), "unseen")
      bad_type = original[1:3,]
      bad_type[["x value"]] = factor(bad_type[["x value"]])
      expect_error(predict(fit, bad_type), "numeric term")
      missing = original[c(7, 2, 9, 1, 8, 3),]
      missing[["x value"]][c(1L, 4L)] = NA_real_
      missing[["off value"]][6L] = NA_real_
      missing_expected = expected_eta[c(7, 2, 9, 1, 8, 3)]
      missing_expected[c(1L, 4L, 6L)] = NA_real_
      missing_actual = predict(fit, missing, type = "link")
      expect_equal(length(missing_actual), nrow(missing))
      expect_identical(is.na(missing_actual), is.na(missing_expected))
      expect_equal(missing_actual, missing_expected, tolerance = 1e-6)
      out = predict(fit, original[rev(seq_len(n)),], conn = conn, row_id = "id", type = "link")
      actual = as.data.frame(out)
      expect_equal(actual[[out[["prediction_column"]]]][match(d[["id"]], actual[["id"]])], expected_eta, tolerance = 1e-6)
      close(out)
      close(fit)
      expect_true(DBI::dbIsValid(conn))
    }
  }

  # Original raw columns survive source projection for arbitrary later cluster
  # requests, sample exports, original positions and response residuals.
  d[["y value"]] = eta + rnorm(n, sd = .5)
  for (retention in c("inference", "all")) {
    fit = dbglm(fml, data = d, weights = "weight value", offset = "off value",
      subset = ~`in sample`, retain = retention, ssc = bare_ssc, strategy = "compress")
    updated = summary(fit, vcov = ~`future cluster`)
    reference = glm.fit(X[d[["in sample"]],,drop = FALSE], d[["y value"]][d[["in sample"]]],
      weights = d[["weight value"]][d[["in sample"]]], offset = d[["off value"]][d[["in sample"]]], family = gaussian())
    expect_equal(as.numeric(vcov(updated)), as.numeric(covariance(X[d[["in sample"]],,drop = FALSE],
      d[["y value"]][d[["in sample"]]], d[["weight value"]][d[["in sample"]]],
      reference[["fitted.values"]], gaussian(), d[["future cluster"]][d[["in sample"]]])), tolerance = 1e-6)
    expect_equal(case.names(fit), which(d[["in sample"]]))
    expect_equal(fit[["offset"]], "off value")
    expect_equal(fit[["weights"]], "weight value")
    expect_equal(fit[["cluster_var"]], NULL)
    expect_identical(fit[["fitted_state"]][["inference"]][["coefficient_ids"]], fit[["model_spec"]][["fitted_ids"]])
    close(fit)
  }

  selector_state = new.env(parent = emptyenv())
  selector_state[["calls"]] = 0L
  reference_level = function() {
    selector_state[["calls"]] = selector_state[["calls"]] + 1L
    1
  }
  selector_fit = dbglm(`y value` ~ `x value` + i(`period value`, ref = reference_level()), data = d)
  expect_identical(selector_state[["calls"]], 1L)
  expect_true(all(is.finite(predict(selector_fit, d))))
  expect_identical(selector_state[["calls"]], 1L)
  close(selector_fit)

  alias_data = d
  alias_data[["DBREG_INTERNAL_X_0001__"]] = d[["x_value"]]
  alias_fit = dbglm(`y value` ~ `x value` + DBREG_INTERNAL_X_0001__, data = alias_data)
  alias_reference = lm.fit(cbind(1, d[["x value"]], d[["x_value"]]), d[["y value"]])
  expect_equal(predict(alias_fit, alias_data, type = "link"), alias_reference[["fitted.values"]], tolerance = 1e-8)
  expect_false(any(tolower(alias_fit[["model_spec"]][["columns"]][["alias"]]) == "dbreg_internal_x_0001__"))
  expect_false(any(tolower(alias_fit[["model_spec"]][["columns"]][["alias"]]) %in% tolower(names(alias_data))))
  expect_identical(anyDuplicated(tolower(alias_fit[["model_spec"]][["columns"]][["alias"]])), 0L)
  broken = alias_fit
  broken[["model_spec"]][["columns"]][["recipe"]] = NULL
  expect_error(predict(broken, alias_data), "Incomplete")
  expect_error(dbreg:::sql_fitted_prediction_design(broken, conn), "Incomplete")
  absent_descriptor = alias_fit
  absent_descriptor[["model_spec"]][["columns"]][["id"]] = NULL
  absent_descriptor[["model_spec"]][["columns"]][["recipe"]] = NULL
  absent_descriptor[["model_spec"]][["source_map"]] = NULL
  expect_error(predict(absent_descriptor, alias_data), "Incomplete")
  expect_error(dbreg:::sql_fitted_prediction_design(absent_descriptor, conn), "Incomplete")
  reordered = alias_fit
  reordered[["model_spec"]][["fitted_ids"]] = rev(reordered[["model_spec"]][["fitted_ids"]])
  expect_error(predict(reordered, alias_data), "completed inference order")
  expect_error(dbreg:::sql_fitted_prediction_design(reordered, conn), "completed inference order")
  repeated = alias_fit
  repeated[["model_spec"]][["fitted_ids"]][2L] = repeated[["model_spec"]][["fitted_ids"]][1L]
  expect_error(predict(repeated, alias_data), "Invalid fitted")
  invalid_sources = alias_fit
  bindings = invalid_sources[["model_spec"]][["source_map"]]
  invalid_sources[["model_spec"]][["source_map"]] = c(bindings, bindings[1L])
  expect_error(predict(invalid_sources, alias_data), "raw source binding")
  close(alias_fit)

  # A shared raw source can be numeric, forced categorical and an i() key.
  # This collinear design checks fitted subset IDs independently of full X.
  mixed = dbglm(`y value` ~ `numeric group` + factor(`numeric group`) +
    i(`numeric group`, `treat value`, ref = 1), data = d, ssc = bare_ssc)
  mixed_X = cbind(1, d[["numeric group"]], as.numeric(d[["numeric group"]] == 2),
    as.numeric(d[["numeric group"]] == 3), as.numeric(d[["numeric group"]] == 2) * d[["treat value"]],
    as.numeric(d[["numeric group"]] == 3) * d[["treat value"]])
  mixed_ref = glm.fit(mixed_X, d[["y value"]], family = gaussian())
  expect_equal(predict(mixed, d, type = "link"), mixed_ref[["linear.predictors"]], tolerance = 1e-8)
  expect_true(length(mixed[["model_spec"]][["fitted_ids"]]) < 1L + nrow(mixed[["model_spec"]][["columns"]]))
  absent = d
  absent[["numeric group"]] = NULL
  bindings = mixed[["model_spec"]][["source_map"]]
  for (name in names(bindings)[bindings == "numeric group"]) absent[[name]] = 1
  expect_error(predict(mixed, absent), "numeric group")
  close(mixed)

  # FE companions preserve complete recovery columns and raw public keys/trends.
  d[["group key"]] = rep(seq_len(20L), length.out = n)
  d[["second key"]] = rep(seq_len(4L), each = 80L)
  d[["time value"]] = rnorm(n)
  d[["y value"]] = .4 * d[["x value"]] + .2 * d[["x_value"]] +
    sin(d[["group key"]]) + .1 * d[["second key"]] + rnorm(n, sd = .5)
  for (fml in list(`y value` ~ `x value` + x_value | `group key`,
    `y value` ~ `x value` + x_value | `group key` + `second key`,
    `y value` ~ `x value` + x_value | `group key`[`time value`])) {
    rhs = if (grepl("time value", paste(deparse(fml), collapse = ""), fixed = TRUE))
      ~ `x value` + x_value + factor(`group key`) * `time value` else
      if (grepl("second key", paste(deparse(fml), collapse = ""), fixed = TRUE))
        ~ `x value` + x_value + factor(`group key`) + factor(`second key`) else
        ~ `x value` + x_value + factor(`group key`)
    reference = lm.wfit(model.matrix(rhs, d), d[["y value"]], d[["weight value"]])
    fit = dbglm(fml, data = d, weights = "weight value", fixef_rm = "none", retain = "all")
    minimal = d[c("x value", "x_value", "group key", "second key", "time value", "id")]
    expect_equal(predict(fit, minimal, type = "link"), reference[["fitted.values"]], tolerance = 1e-6)
    partial = minimal[c(7, 2, 9, 1, 8, 3),]
    partial[["x value"]][c(1L, 4L)] = NA_real_
    partial[["group key"]][3L] = 999L
    expected_partial = as.numeric(reference[["fitted.values"]][c(7, 2, 9, 1, 8, 3)])
    expected_partial[c(1L, 3L, 4L)] = NA_real_
    actual_partial = predict(fit, partial, type = "link")
    expect_equal(length(actual_partial), nrow(partial))
    expect_identical(is.na(actual_partial), is.na(expected_partial))
    expect_equal(actual_partial, expected_partial, tolerance = 1e-6)
    expect_identical(fit[["fitted_state"]][["specification"]][["dictionaries"]][[1L]][["columns"]], "group key")
    out = predict(fit, minimal, conn = conn, row_id = "id", type = "link")
    actual = as.data.frame(out)
    expect_equal(actual[[out[["prediction_column"]]]][match(d[["id"]], actual[["id"]])], reference[["fitted.values"]], tolerance = 1e-6)
    close(out)
    close(fit)
  }
  only_fe = dbglm(`y value` ~ 1 | `group key`, data = d, weights = "weight value", retain = "all", fixef_rm = "none")
  means = ave(seq_len(n), d[["group key"]], FUN = function(rows) weighted.mean(d[["y value"]][rows], d[["weight value"]][rows]))
  expect_equal(predict(only_fe, d["group key"], type = "link"), as.numeric(means), tolerance = 1e-6)
  expect_identical(only_fe[["model_spec"]][["fitted_ids"]], character())
  expect_identical(attr(vcov(only_fe), "dbreg_coefficient_ids"), character())
  expect_identical(only_fe[["fitted_state"]][["inference"]][["coefficient_ids"]], character())
  close(only_fe)

  copy_data = d
  copy_data[["x copy"]] = d[["x value"]]
  fit = dbglm(`y value` ~ `x value` + `x copy` | `group key`, data = copy_data,
    weights = "weight value", fixef_rm = "none", retain = "prediction")
  reference = lm.wfit(model.matrix(~ `x value` + factor(`group key`), d), d[["y value"]], d[["weight value"]])
  expect_equal(predict(fit, copy_data, type = "link"), reference[["fitted.values"]], tolerance = 1e-6)
  expect_equal(nrow(fit[["model_spec"]][["columns"]]), 2)
  expect_equal(length(fit[["model_spec"]][["fitted_ids"]]), 1)
  expect_equal(length(fit[["fitted_state"]][["specification"]][["design_names"]]), 2)
  close(fit)
  old = options(dbreg.prediction_max_entries = 32)
  on.exit(options(old), add = TRUE)
  fit = dbglm(`y value` ~ `x value` + x_value | `group key`, data = d,
    weights = "weight value", fixef_rm = "none", retain = "prediction")
  options(old)
  expect_identical(fit[["fitted_state"]][["prediction"]][["storage"]], "database")
  reference = lm.wfit(model.matrix(~ `x value` + x_value + factor(`group key`), d), d[["y value"]], d[["weight value"]])
  expect_equal(predict(fit, d, type = "link"), reference[["fitted.values"]], tolerance = 1e-6)
  partial = d[c(7, 2, 9, 1, 8, 3),]
  partial[["x value"]][c(1L, 4L)] = NA_real_
  partial[["group key"]][3L] = 999L
  expected_partial = as.numeric(reference[["fitted.values"]][c(7, 2, 9, 1, 8, 3)])
  expected_partial[c(1L, 3L, 4L)] = NA_real_
  actual_partial = predict(fit, partial, type = "link")
  expect_equal(length(actual_partial), nrow(partial))
  expect_identical(is.na(actual_partial), is.na(expected_partial))
  expect_equal(actual_partial, expected_partial, tolerance = 1e-6)
  expect_error(predict(fit, d, conn = conn, row_id = "id", type = "link"), "owning session")
  out = predict(fit, d, row_id = "id", type = "link")
  actual = as.data.frame(out)
  expect_equal(actual[[out[["prediction_column"]]]][match(d[["id"]], actual[["id"]])], reference[["fitted.values"]], tolerance = 1e-6)
  close(out)
  close(fit)
  expect_error(predict(fit, d), "released|closed")

  DBI::dbWriteTable(conn, "glm_source_rows", d)
  DBI::dbExecute(conn, 'ALTER TABLE glm_source_rows RENAME COLUMN id TO "row key"')
  fit = dbglm(`y value` ~ `x value` + x_value, conn = conn, table = "glm_source_rows",
    subset = ~`in sample`, row_id = "row key", sample_table = "glm_raw_keys", retain = "inference")
  expect_identical(DBI::dbListFields(conn, "glm_raw_keys"), "row key")
  keys = DBI::dbReadTable(conn, "glm_raw_keys", check.names = FALSE)
  expect_identical(names(keys), "row key")
  expect_equal(sort(as.numeric(keys[["row key"]])), which(d[["in sample"]]))
  close(fit)
  expect_true(DBI::dbExistsTable(conn, "glm_source_rows"))
  expect_true(DBI::dbExistsTable(conn, "glm_raw_keys"))
})
