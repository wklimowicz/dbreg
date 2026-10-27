# Labels identify fitted sources/bases; literal source/level text is opaque.
local({
  set.seed(11061)
  d = data.frame(a = rnorm(180), b = rnorm(180), a_x_b = rnorm(180),
    f_x_cat = sample(c("a", "b_x_c", "c"), 180, TRUE), z_x = rep(1:3, 60))
  d[["y"]] = 1 + d[["a"]] - 0.7 * d[["a_x_b"]] + 0.4 * d[["a"]] * d[["b"]] + rnorm(180)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "structured_label_training", d)
  direct = dbreg::sql_model_matrix(y ~ a_x_b + a:b, conn, "structured_label_training")
  expect_identical(direct[["col_names"]], c("a_x_b", "a_x_b"))
  expect_identical(direct[["aliases"]][[1L]], "a_x_b")
  expect_true(direct[["aliases"]][[2L]] != "a_x_b")
  expect_identical(anyDuplicated(tolower(direct[["aliases"]])), 0L)
  expect_identical(direct[["coef_names"]], c("a_x_b", "a:b"))
  expect_identical(direct[["columns"]][["label"]], direct[["coef_names"]])
  formulas = list(y ~ a_x_b, y ~ a:b + a_x_b, y ~ a_x_b + f_x_cat,
    y ~ a:f_x_cat + a_x_b, y ~ a * factor(z_x), y ~ I(a_x_b^2) + f_x_cat,
    y ~ 0 + factor(z_x), y ~ 0 + f_x_cat)
  for (formula in formulas) {
    reference = stats::lm(formula, data = d)
    X = stats::model.matrix(reference)
    beta = as.numeric(qr.solve(X, d[["y"]]))
    V = solve(crossprod(X)) * sum((d[["y"]] - X %*% beta)^2) / (nrow(d) - ncol(X))
    rows = d[c(11, 4, 20, 1, 31), ]
    rows[["id"]] = c(15L, 2L, 99L, 1L, 8L)
    expected = as.numeric(stats::predict(reference, rows))
    for (strategy in c("moments", "compress")) {
      model = dbreg::dbreg(formula, data = d, conn = conn, strategy = strategy)
      expect_identical(names(coef(model)), names(stats::coef(reference)))
      expect_equal(as.numeric(coef(model)), beta, tolerance = 1e-9)
      expect_equal(as.numeric(vcov(model)), as.numeric(V), tolerance = 1e-9)
      expect_identical(rownames(vcov(model)), names(coef(model)))
      expect_identical(generics::tidy(model)[["term"]], names(coef(model)))
      expect_equal(predict(model, rows), expected, tolerance = 1e-9)
      expect_equal(as.numeric(as.matrix(predict(model, rows, interval = "confidence"))),
        as.numeric(stats::predict(reference, rows, interval = "confidence")), tolerance = 1e-9)
      extra = rows
      extra[["a_x_b_1"]] = 999
      extra[["dbreg_internal_x_0001__"]] = -999
      expect_equal(predict(model, extra), expected, tolerance = 1e-9)
      out = predict(model, extra[c(5, 1, 3, 2, 4), ], conn = conn, row_id = "id")
      actual = as.data.frame(out)
      expect_equal(actual[[out[["prediction_column"]]]][match(rows[["id"]], actual[["id"]])], expected, tolerance = 1e-9)
      close(out)
      close(model)
    }
  }
  for (formula in list(y ~ a_x_b + f_x_cat, y ~ a_x_b + a:b)) {
    reference = stats::lm(formula, data = d)
    model = dbreg::dbglm(formula, data = d, family = stats::gaussian(), retain = "prediction")
    expect_identical(names(coef(model)), names(coef(reference)))
    expect_equal(as.numeric(coef(model)), as.numeric(coef(reference)), tolerance = 1e-8)
    expect_equal(predict(model, d, type = "response"), as.numeric(predict(reference, d)), tolerance = 1e-8)
    close(model)
  }
  # Typed i() keeps its existing convention,including literal level text.
  X = cbind(1, as.numeric(d[["f_x_cat"]] == "b_x_c"), as.numeric(d[["f_x_cat"]] == "c"), d[["a_x_b"]])
  beta = as.numeric(qr.solve(X, d[["y"]]))
  for (strategy in c("moments", "compress")) {
    model = dbreg::dbreg(y ~ i(f_x_cat, ref = "a") + a_x_b, data = d, conn = conn, strategy = strategy)
    expect_identical(names(coef(model)), c("(Intercept)", "f_x_cat::b_x_c", "f_x_cat::c", "a_x_b"))
    expect_equal(as.numeric(coef(model)), beta, tolerance = 1e-9)
    expect_equal(predict(model, d), as.numeric(X %*% beta), tolerance = 1e-9)
    close(model)
    one = d
    one[["only_x_value"]] = 1
    model = dbreg::dbreg(y ~ 0 + factor(only_x_value), data = one, conn = conn, strategy = strategy)
    expect_identical(names(coef(model)), "factor(only_x_value)_constant")
    expect_equal(as.numeric(coef(model)), mean(one[["y"]]), tolerance = 1e-9)
    expect_equal(predict(model, one[1:3, ]), rep(mean(one[["y"]]), 3), tolerance = 1e-9)
    close(model)
  }
})

# Cross-block IV targets are distinct; public names retain the actual sources.
local({
  set.seed(11062)
  d = data.frame(a = rnorm(180), b_x_c = rnorm(180), a_x_b = rnorm(180), c = rnorm(180))
  d[["endog"]] = 0.7 * d[["a_x_b"]] * d[["c"]] + rnorm(180)
  d[["y"]] = 1 + 0.4 * d[["a"]] * d[["b_x_c"]] + 1.2 * d[["endog"]] + rnorm(180)
  X = cbind(1, d[["a"]] * d[["b_x_c"]], d[["endog"]])
  Z = cbind(1, X[, 2L], d[["a_x_b"]] * d[["c"]])
  bread = solve(crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, X)))
  beta = as.numeric(bread %*% crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, d[["y"]])))
  V = bread * sum((d[["y"]] - X %*% beta)^2) / 177
  for (strategy in c("moments", "compress")) {
    model = dbreg::dbivreg(y ~ a:b_x_c | endog ~ a_x_b:c, data = d, strategy = strategy)
    expect_identical(names(coef(model)), c("(Intercept)", "a:b_x_c", "endog"))
    expect_equal(as.numeric(coef(model)), beta, tolerance = 1e-9)
    expect_equal(as.numeric(vcov(model)), as.numeric(V), tolerance = 1e-9)
    expect_equal(attr(vcov(model), "df.K"), 3)
    expect_equal(predict(model, d), as.numeric(X %*% beta), tolerance = 1e-9)
    close(model)
  }
})

# A nonstandard factor's dummy basename cannot rename a raw numeric source.
local({
  set.seed(11067)
  d = data.frame(grp_labelb = rnorm(180), z = rnorm(180))
  d[["grp label"]] = rep(c("a", "b"), 90)
  d[["endog"]] = 0.8 * d[["z"]] + rnorm(180)
  d[["y"]] = 1 + 0.4 * d[["grp_labelb"]] + (d[["grp label"]] == "b") + 1.2 * d[["endog"]] + rnorm(180)
  X = cbind(1, as.numeric(d[["grp label"]] == "b"), d[["grp_labelb"]], d[["endog"]])
  Z = cbind(X[, 1:3], d[["z"]])
  bread = solve(crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, X)))
  beta = as.numeric(bread %*% crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, d[["y"]])))
  for (strategy in c("moments", "compress")) {
    model = dbreg::dbivreg(y ~ `grp label` + grp_labelb | endog ~ z, data = d, strategy = strategy)
    expect_identical(names(coef(model)), c("(Intercept)", "`grp label`b", "grp_labelb", "endog"))
    expect_equal(as.numeric(coef(model)), beta, tolerance = 1e-9)
    expect_equal(predict(model, d), as.numeric(X %*% beta), tolerance = 1e-9)
    close(model)
  }
})

# Learned effects and retained axes preserve source labels and typed FE coding.
local({
  set.seed(11063)
  d = data.frame(id = 1:240, x_x_raw = rnorm(240), f_x_cat = sample(c("a", "b_x_c", "c"), 240, TRUE),
    g_x_id = rep(1:8, each = 30), z_x = rnorm(240), cl = rep(1:12, 20))
  d[["y"]] = 1 + d[["x_x_raw"]] + sin(d[["g_x_id"]]) * d[["z_x"]] + rnorm(240)
  reference = stats::lm(y ~ x_x_raw + factor(f_x_cat) + factor(g_x_id) + factor(g_x_id):z_x, data = d)
  for (strategy in c("compress", "demean")) {
    model = dbreg::dbreg(y ~ x_x_raw + factor(f_x_cat) | g_x_id[z_x], data = d,
      strategy = strategy, retain = "all", vcov = ~cl)
    expected_names = c("x_x_raw", "factor(f_x_cat)b_x_c", "factor(f_x_cat)c")
    expect_identical(names(coef(model)), expected_names)
    expect_equal(as.numeric(coef(model)), as.numeric(coef(reference)[2:4]), tolerance = 1e-8)
    expect_equal(predict(model, d), as.numeric(stats::predict(reference, d)), tolerance = 1e-8)
    expect_equal(predict(model), as.numeric(stats::predict(reference, d)), tolerance = 1e-8)
    expect_identical(rownames(vcov(model)), expected_names)
    changed = summary(model, vcov = "hc1")
    expect_identical(names(coef(changed)), expected_names)
    expect_identical(generics::tidy(changed)[["term"]], expected_names)
    close(model)
  }
})

# Mundlak retains its current means annotation while preserving literal sources.
local({
  set.seed(11064)
  d = data.frame(x_x_raw = rnorm(240), f_x_cat = sample(c("a", "b_x_c", "c"), 240, TRUE), g_x_id = rep(1:8, each = 30))
  d[["y"]] = 1 + d[["x_x_raw"]] + sin(d[["g_x_id"]]) + rnorm(240)
  source = stats::model.matrix(~x_x_raw + f_x_cat, d)
  means = apply(source[, -1, drop = FALSE], 2L, function(x) ave(x, d[["g_x_id"]], FUN = mean))
  X = cbind(source, means)
  beta = as.numeric(qr.solve(X, d[["y"]]))
  V = solve(crossprod(X)) * sum((d[["y"]] - X %*% beta)^2) / (nrow(d) - ncol(X))
  model = dbreg::dbreg(y ~ x_x_raw + f_x_cat | g_x_id, data = d, strategy = "mundlak")
  labels = c("x_x_raw", "f_x_catb_x_c", "f_x_catc")
  expect_identical(names(coef(model)), labels)
  expect_identical(names(coef(model, fe = TRUE)), c("(Intercept)", labels, paste0(labels, "_bar_g_x_id")))
  expect_equal(as.numeric(coef(model, fe = TRUE)), beta, tolerance = 1e-8)
  expect_equal(as.numeric(vcov(model, fe = TRUE)), as.numeric(V), tolerance = 1e-8)
  expect_equal(predict(model, d), as.numeric(X %*% beta), tolerance = 1e-8)
  close(model)
})

# Bins outer source mapping re-renders typed descriptors,including FE controls.
local({
  set.seed(11065)
  d = data.frame(x = stats::runif(240))
  d[["control_x label"]] = sample(c("a", "b_x_c", "c"), 240, TRUE)
  d[["group_x label"]] = rep(1:8, each = 30)
  d[["y"]] = 1 + 2 * d[["x"]] + sin(d[["group_x label"]]) + rnorm(240)
  breaks = seq(0, 1, length.out = 5)
  d[["b"]] = cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE)
  for (with_fe in c(FALSE, TRUE)) {
    formula = if (with_fe) y ~ x + `control_x label` | `group_x label` else y ~ x + `control_x label`
    reference_formula = if (with_fe) y ~ b + `control_x label` + factor(`group_x label`) else y ~ b + `control_x label`
    reference = stats::lm(reference_formula, d)
    for (strategy in c("auto", "compress")) {
      output = dbreg::dbbinsreg(formula, data = d, nbins = 4, binspos = breaks,
        points = c(0, 0), strategy = strategy, vcov = "iid", plot = FALSE)
      model = output[["model"]]
      expect_identical(model[["inference_result"]][["kind"]], "iid")
      expect_equal(as.numeric(coef(model, fe = TRUE)), as.numeric(coef(reference)), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(model, fe = TRUE)), as.numeric(vcov(reference)), tolerance = 1e-8)
      expect_identical(names(coef(model, fe = TRUE))[5:6], c("`control_x label`b_x_c", "`control_x label`c"))
      expect_identical(rownames(vcov(model, fe = TRUE)), names(coef(model, fe = TRUE)))
      expect_identical(model[["inference_result"]][["coefficient_names"]], rownames(model[["coeftable"]]))
      expect_identical(names(coef(model, collin = TRUE)), names(coef(model)))
      close(model)
    }
  }
})

# Actual dense/sparse producer coding,including no-intercept promotion and FE
# nuisance interactions,determines recipe order without inspecting label text.
local({
  set.seed(11068)
  for (large_basis in c(FALSE, TRUE)) {
    n = if (large_basis) 240L else 36L
    f = if (large_basis) sprintf("f_x_%02d", 1:30) else c("a", "b_x_c", "c")
    g = if (large_basis) sprintf("g_x_%02d", 1:40) else c("a", "b_x_c")
    d = data.frame(a_x_b = rnorm(n), f_x_cat = factor(rep(f, length.out = n), levels = f),
      g_x_cat = factor(rep(g, length.out = n), levels = g))
    sparse_seen = FALSE
    formulas = list(~0 + f_x_cat + g_x_cat, ~0 + a_x_b + f_x_cat + g_x_cat,
      ~0 + g_x_cat + g_x_cat:a_x_b, ~0 + f_x_cat + g_x_cat + g_x_cat:a_x_b,
      ~f_x_cat + g_x_cat + g_x_cat:a_x_b, ~0 + a_x_b:f_x_cat, ~0 + g_x_cat:f_x_cat)
    for (formula in formulas) {
      X = dbreg:::fitted_model_matrix(formula, d)
      if (inherits(X, "Matrix")) sparse_seen = TRUE
      recipes = dbreg:::model_matrix_recipes(X, formula, d)
      labels = vapply(recipes, dbreg:::linear_recipe_label, character(1))
      contrasts = list(f_x_cat = stats::contr.treatment(f), g_x_cat = stats::contr.treatment(g))
      reference = stats::model.matrix(formula, d, contrasts.arg = contrasts[intersect(names(contrasts), all.vars(formula))])
      expect_identical(labels, colnames(reference))
      values = do.call(cbind, lapply(recipes, dbreg:::linear_recipe_value, data = d, i_specs = list()))
      expect_equal(as.numeric(values), as.numeric(reference), tolerance = 1e-14)
      expect_equal(as.numeric(X), as.numeric(reference), tolerance = 1e-14)
    }
    if (large_basis) expect_true(sparse_seen)
  }
})

# Retained FE-only fits carry an explicit empty slope/name/ID axis.
local({
  set.seed(1827)
  n = 240L
  d = data.frame(id = seq_len(n), g = rep(1:8, length.out = n),
    h = sample.int(4L, n, TRUE), j = sample.int(3L, n, TRUE))
  d[["y"]] = sin(d[["g"]]) + cos(d[["h"]]) + stats::rnorm(n)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  for (third in c(FALSE, TRUE)) {
    formula = if (third) y ~ 1 | g + h + j else y ~ 1 | g + h
    reference = fixest::feols(formula, d, fixef.rm = "none")
    X = stats::model.matrix(if (third) ~factor(g) + factor(h) + factor(j) else ~factor(g) + factor(h), d)
    beta = as.numeric(qr.solve(X, d[["y"]]))
    K = ncol(X)
    V = solve(crossprod(X)) * sum((d[["y"]] - X %*% beta)^2) / (n - K)
    for (strategy in c("compress", "demean", "auto")) {
      model = dbreg::dbreg(formula, data = d, conn = conn, strategy = strategy,
        retain = "all", vcov = "iid", ssc = list(K.fixef = "full"), fixef_rm = "none")
      expect_identical(length(stats::coef(model)), 0L)
      expect_identical(dim(stats::vcov(model)), c(0L, 0L))
      expect_identical(model[["linear_design"]][["slope_ids"]], character())
      expect_identical(model[["fitted_state"]][["inference"]][["coefficient_names"]], character())
      expect_identical(model[["fitted_state"]][["inference"]][["coefficient_ids"]], character())
      expect_equal(attr(stats::vcov(model), "df.K"), K)
      expect_equal(stats::df.residual(model), n - K)
      expect_equal(stats::predict(model, d), as.numeric(X %*% beta), tolerance = 1e-7)
      expect_equal(stats::predict(model, d), as.numeric(stats::predict(reference, d)), tolerance = 1e-7)
      if (identical(model[["strategy"]], "compress")) {
        expect_equal(as.numeric(stats::coef(model, fe = TRUE)), beta, tolerance = 1e-7)
        expect_equal(as.numeric(stats::vcov(model, fe = TRUE)), as.numeric(V), tolerance = 1e-7)
      }
      rows = d[c(11, 4, 20, 1, 31), ]
      output = stats::predict(model, newdata = rows, conn = conn, row_id = "id")
      actual = as.data.frame(output)
      expected = as.numeric(X %*% beta)[match(actual[["id"]], d[["id"]])]
      expect_equal(actual[[output[["prediction_column"]]]], expected, tolerance = 1e-7)
      close(output)
      close(model)
    }
  }
})
