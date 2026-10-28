# Independent source/dictionary role witnesses. Labels never define the oracle X.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  empty_inputs = list2env(list(fml = Formula::Formula(y ~ 1), factor_levels = list(),
    factor_sources = list(), linear_source_map = c(y = "y"), source_classes = character(),
    fe = NULL, fe_slopes = list(), source_map = character()), parent = emptyenv())
  empty_layout = dbreg:::complete_linear_source_contract(list(), empty_inputs)
  for (field in c("fe_sources", "trend_sources", "factor_roles")) expect_identical(unname(empty_layout[[field]]), character())
  index = seq_len(240)
  d = data.frame(id = index, y = 0, f = rep(c(1, 1, 2, 2), 60),
    literal_source = rep(c("a", "b"), 120), g = rep(seq_len(40), each = 6),
    w = 1 + (index %% 5) / 4)
  d[["y"]] = 1 + .7 * (d[["f"]] == 2) + 1.3 * (d[["literal_source"]] == "b") +
    .2 * sin(d[["g"]]) + .05 * sin(index * .71)
  check_ols = function(strategy, collision = TRUE, fe = FALSE, database = FALSE, composite = FALSE, retain = if (fe) "prediction" else "none") {
    old = options(dbreg.prediction_max_entries = if (database) 64 else 100000)
    on.exit(options(old))
    data = d
    if (composite) {
      data[["g"]] = rep(seq_len(10), each = 24)
      data[["h"]] = rep(rep(seq_len(4), each = 6), 10)
    }
    literal = if (collision) "factor(f)" else "literal_source"
    names(data)[names(data) == "literal_source"] = literal
    indicator = function(rows) cbind(as.numeric(rows[[literal]] == "b"), as.numeric(rows[["f"]] == 2))
    matrix = function(rows) cbind(indicator(rows), if (fe) outer(if (composite) (rows[["g"]] - 1L) * 4L + rows[["h"]] else rows[["g"]], seq_len(40), `==`) * 1 else 1)
    X = matrix(data)
    reference = stats::lm.wfit(X, data[["y"]], data[["w"]])
    expect_equal(reference[["rank"]], ncol(X))
    V = solve(crossprod(X, data[["w"]] * X)) *
      sum(data[["w"]] * reference[["residuals"]]^2) / (nrow(X) - ncol(X))
    first = paste(deparse(as.name(literal), backtick = TRUE), collapse = "")
    formula = stats::as.formula(paste("y ~", first, "+ factor(f)", if (composite) "| g ^ h" else if (fe) "| g" else ""))
    model = dbreg(formula, data = data, conn = conn, weights = "w", strategy = strategy,
      retain = retain, vcov = "iid", ssc = list(K.fixef = "full"))
    on.exit(close(model), add = TRUE, after = FALSE)
    # The independent X puts its intercept last; public coefficients put it first.
    selected = if (fe) 1:2 else c(3L, 1L, 2L)
    expect_equal(as.numeric(coef(model)), as.numeric(reference[["coefficients"]][selected]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(model)), as.numeric(V[selected, selected, drop = FALSE]), tolerance = 1e-7)
    layout = model[["linear_design"]]
    expect_identical(typeof(layout[["fe_sources"]]), "character")
    expect_identical(unname(layout[["trend_sources"]]), character())
    if (!fe) expect_identical(unname(layout[["fe_sources"]]), character())
    if (fe) {
      expect_true(length(layout[["fe_sources"]]) > 0L)
      missing_fe = model
      missing_fe[["linear_design"]][["fe_sources"]] = NULL
      expect_error(predict(missing_fe, data), "descriptor|source")
      expect_error(predict(missing_fe, data, conn = conn, row_id = "id"), "descriptor|source")
      missing_fe[["linear_design"]][["fe_sources"]] = layout[["fe_sources"]][-1L]
      expect_error(predict(missing_fe, data), "descriptor|source")
      expect_error(predict(missing_fe, data, conn = conn, row_id = "id"), "descriptor|source")
    }
    expect_false(anyDuplicated(names(model[["factor_levels"]])) > 0L)
    expect_true(all(c("ordinary", "forced") %in% layout[["factor_roles"]]))
    rows = data[c(17, 2, 201, 55, 120, 9, 37, 89, 160, 230), ]
    rows[[literal]] = ifelse(rows[[literal]] == "a", "b", "a")
    rows[["f"]] = 3 - rows[["f"]]
    rows[["y"]] = -999
    for (alias in setdiff(names(layout[["source_map"]]), unname(layout[["source_map"]]))) rows[[alias]] = "fake"
    for (alias in names(model[["composite_fes"]])) rows[[alias]] = "fake-derived-rank"
    expected = as.numeric(matrix(rows) %*% reference[["coefficients"]])
    expect_equal(as.numeric(predict(model, rows)), expected, tolerance = 1e-7)
    if (fe && retain == "none") {
      expect_error(predict(model, rows, conn = conn, row_id = "id"), "retain")
    } else {
      out = predict(model, rows, conn = conn, row_id = "id")
      table = as.data.frame(out)
      actual = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
      close(out)
      expect_equal(as.numeric(actual), expected, tolerance = 1e-7)
    }
    if (fe && retain != "none") expect_identical(identical(model[["fitted_state"]][["prediction"]][["storage"]], "database"), database)
    bad = rows
    bad[[literal]][1L] = "unseen"
    expect_error(predict(model, bad), "unseen|level")
    bad = rows
    bad[["f"]][1L] = 3
    expect_error(predict(model, bad), "unseen|level")
    broken = model
    broken[["linear_design"]][["source_map"]] = NULL
    expect_error(predict(broken, rows), "descriptor|bound|source")
    expect_error(predict(broken, rows, conn = conn, row_id = "id"), "descriptor|bound|source")
    if (!fe) {
      built = dbreg:::native_linear_matrix(model, rows)
      expect_null(attr(built, "dbreg_column_ids", exact = TRUE))
    }
  }
  check_ols("moments", collision = FALSE)
  check_ols("moments")
  check_ols("compress")
  check_ols("demean", fe = TRUE)
  check_ols("demean", fe = TRUE, database = TRUE)
  check_ols("compress", fe = TRUE)
  check_ols("demean", fe = TRUE, composite = TRUE)
  check_ols("compress", fe = TRUE, composite = TRUE, retain = "none")

  # A raw FE key and a forced slope dictionary can have identical rendered text.
  # Independent group-dummy X proves both key and slope uploads bind raw roles.
  data = d
  names(data)[names(data) == "literal_source"] = "factor(f)"
  full_X = function(rows) cbind(as.numeric(rows[["f"]] == 2),
    as.numeric(rows[["factor(f)"]] == "a"), as.numeric(rows[["factor(f)"]] == "b"))
  X = full_X(data)
  reference = stats::lm.wfit(X, data[["y"]], data[["w"]])
  V = solve(crossprod(X, data[["w"]] * X)) * sum(data[["w"]] * reference[["residuals"]]^2) / (nrow(X) - ncol(X))
  fe_model = dbreg(y ~ factor(f) | `factor(f)`, data = data, conn = conn, weights = "w",
    strategy = "demean", retain = "prediction", ssc = list(K.fixef = "full"))
  expect_equal(as.numeric(coef(fe_model)), as.numeric(reference[["coefficients"]][1L]), tolerance = 1e-7)
  expect_equal(as.numeric(vcov(fe_model)), as.numeric(V[1L, 1L]), tolerance = 1e-7)
  rows = data[c(17, 2, 55, 9, 37, 89), ]
  rows[["factor(f)"]] = ifelse(rows[["factor(f)"]] == "a", "b", "a")
  rows[["f"]] = 3 - rows[["f"]]
  expected = as.numeric(full_X(rows) %*% reference[["coefficients"]])
  expect_equal(as.numeric(predict(fe_model, rows)), expected, tolerance = 1e-7)
  prepared = dbreg:::prediction_dataframe_source(fe_model, rows, c("id", "f", "factor(f)"))
  forced = names(fe_model[["linear_design"]][["factor_roles"]])[fe_model[["linear_design"]][["factor_roles"]] == "forced"]
  expect_false(identical(prepared[["labels"]][[forced]], prepared[["key_labels"]][["factor(f)"]]))
  out = predict(fe_model, rows, conn = conn, row_id = "id")
  table = as.data.frame(out)
  actual = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
  close(out)
  expect_equal(as.numeric(actual), expected, tolerance = 1e-7)

  close(fe_model)

  # A renamed source used numerically and in i() stays numeric during compression.
  local({
    values = rep(seq_len(4), 60)
    typed = data.frame(id = index, y = .4 * values + .6 * (values == 2) + sin(index * .31), w = d[["w"]])
    typed[["elapsed value"]] = values
    means = vapply(seq_len(4), function(value) {
      rows = values == value
      sum(typed[["w"]][rows] * typed[["y"]][rows]) / sum(typed[["w"]][rows])
    }, numeric(1))
    model = dbreg(y ~ `elapsed value` + i(`elapsed value`, ref = 1), data = typed,
      conn = conn, weights = "w", strategy = "compress")
    on.exit(close(model))
    rows = typed[c(4, 2, 1, 3), ]
    expected = means[rows[["elapsed value"]]]
    expect_equal(as.numeric(predict(model, rows)), expected, tolerance = 1e-7)
    spec = model[["i_specs"]][[1L]]
    expect_identical(spec[["source_type"]], "numeric")
    expect_false(identical(spec[["first"]], "elapsed value"))
    expect_identical(model[["linear_design"]][["source_map"]][[spec[["first"]]]], "elapsed value")
    out = predict(model, rows, conn = conn, row_id = "id")
    table = as.data.frame(out)
    actual = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
    close(out)
    expect_equal(as.numeric(actual), expected, tolerance = 1e-7)
  })

  # Completed IV prediction owns structural X dictionaries; Z stays estimation-only.
  data[["z"]] = sin(index * .37) + cos(index * .19)
  data[["x"]] = .9 * data[["z"]] + .1 * cos(index * .43)
  X = cbind(1, as.numeric(data[["factor(f)"]] == "b"), as.numeric(data[["f"]] == 2), data[["x"]])
  Z = cbind(X[, 1:3], data[["z"]])
  data[["y"]] = as.numeric(X %*% c(1, 1.3, .7, 1.8)) + .03 * sin(index * .79)
  projected = Z %*% solve(crossprod(Z), crossprod(Z, X))
  beta = solve(crossprod(X, projected), crossprod(projected, data[["y"]]))
  V = solve(crossprod(X, projected)) * sum((data[["y"]] - X %*% beta)^2) / (nrow(X) - ncol(X))
  # Forced factor() exogenous terms remain unsupported by the IV parser.
  expect_error(dbivreg(y ~ `factor(f)` + factor(f) | x ~ z, data = data, conn = conn), "Transformed|function")
  observed = data
  observed[["f"]] = factor(observed[["f"]], levels = c(1, 2))
  for (strategy in c("moments", "compress")) {
    model = dbivreg(y ~ `factor(f)` + f | x ~ z, data = observed, conn = conn,
      strategy = strategy, vcov = "iid")
    expect_equal(as.numeric(coef(model)), as.numeric(beta), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(model)), as.numeric(V), tolerance = 1e-7)
    rows = data[c(17, 2, 55, 9, 37, 89), ]
    rows[["z"]] = NULL
    rows[["x"]] = rows[["x"]] + .4
    expected = as.numeric(cbind(1, as.numeric(rows[["factor(f)"]] == "b"), as.numeric(rows[["f"]] == 2), rows[["x"]]) %*% beta)
    expect_equal(as.numeric(predict(model, rows)), expected, tolerance = 1e-7)
    expect_error(predict(model, rows, conn = conn, row_id = "id"), "IV")
    close(model)
  }
})
