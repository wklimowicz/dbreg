# A sole categorical level has the package's declared X = 1 basis. The
# independent oracle supplies that column explicitly: lm(cat) rejects it.
local({
  set.seed(27271)
  n = 96L
  d = data.frame(x = rep(seq(-1, 1, length.out = 8L), 12L),
    w = rep(c(.4, 1, 2, 3), 24L), cl = rep(seq_len(12L), each = 8L),
    h = rep(letters[1:3], 32L), g = rep("sole FE", n),
    off = rep(c(-.3, .2), 48L))
  d[["y"]] = 1.4 + .7 * d[["x"]] + rnorm(n) + d[["off"]]
  sandwich = function(X, y, w, cluster, kind) {
    fit = lm.wfit(X, y, w)
    keep = fit[["qr"]][["pivot"]][seq_len(fit[["rank"]])]
    X = X[, keep, drop = FALSE]
    bread = solve(crossprod(X, w * X))
    e = y - as.numeric(X %*% fit[["coefficients"]][keep])
    N = length(y)
    K = ncol(X)
    V = if (kind == "iid") bread * sum(w * e^2) / (N - K) else {
      scores = X * (w * e)
      if (kind == "cluster") scores = rowsum(scores, cluster, reorder = FALSE)
      meat = crossprod(scores)
      adjustment = if (kind == "hc1") N / (N - K) else
        nrow(scores) / (nrow(scores) - 1) * (N - 1) / (N - K)
      adjustment * bread %*% meat %*% bread
    }
    list(beta = fit[["coefficients"]][keep], V = V,
      prediction = y - e, K = K, keep = keep)
  }
  values = list(character = rep("single", n), literal = rep("_constant", n),
    factor = factor(rep("single", n)), true = rep(TRUE, n), false = rep(FALSE, n))
  for (value in values) for (intercept in c(TRUE, FALSE)) {
    d[["cat"]] = value
    fml = if (intercept) y ~ cat else y ~ 0 + cat
    for (kind in c("iid", "hc1", "cluster")) {
      fit = dbreg(fml, data = d, weights = "w", offset = "off",
        vcov = if (kind == "cluster") ~cl else kind, strategy = "compress")
      oracle = sandwich(matrix(1, n, 1), d[["y"]] - d[["off"]], d[["w"]], d[["cl"]], kind)
      expect_equal(unname(coef(fit)), unname(oracle[["beta"]]), tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(oracle[["V"]]), tolerance = 1e-9)
      expect_equal(fit[["inference_result"]][["n_params_vcov"]], 1)
      expect_equal(df.residual(fit), n - 1)
      expect_equal(names(coef(fit)), if (intercept) "(Intercept)" else "cat_constant")
      expect_equal(predict(fit, d), oracle[["prediction"]] + d[["off"]], tolerance = 1e-9)
      raw = dbreg(fml, data = d, weights = "w", offset = "off",
        strategy = "compress", data_only = TRUE)
      expect_identical(class(fit[["data"]][["cat"]]), class(raw[["cat"]]))
      expect_equal(unique(as.character(fit[["data"]][["cat"]])), unique(as.character(value)))
      nd = d[1:3, , drop = FALSE]
      nd[["cat"]] = if (is.logical(value)) !value[1:3] else rep("unseen", 3L)
      expect_error(predict(fit, nd), "unseen")
      close(fit)
    }
  }

  # The original raw source and numeric source with a constant-looking name
  # have separate fitted coordinates. Verify complete covariance against QR.
  d[["cat"]] = rep("single", n)
  d[["cat_constant"]] = d[["x"]]^2
  for (intercept in c(TRUE, FALSE)) {
    fml = if (intercept) y ~ cat + x + h + cat_constant else y ~ 0 + cat + x + h + cat_constant
    # Recoding the categorical constant to numeric 1 is the declared basis;
    # all other categorical treatment/full-indicator coding remains R's.
    numeric_data = d
    numeric_data[["cat"]] = 1
    X = model.matrix(fml, numeric_data)
    for (kind in c("iid", "hc1", "cluster")) {
      fit = dbreg(fml, data = d, weights = "w", strategy = "compress",
        vcov = if (kind == "cluster") ~cl else kind)
      oracle = sandwich(X, d[["y"]], d[["w"]], d[["cl"]], kind)
      expect_equal(unname(coef(fit)), unname(oracle[["beta"]]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit)), as.numeric(oracle[["V"]]), tolerance = 1e-8)
      expect_equal(predict(fit, d), oracle[["prediction"]], tolerance = 1e-8)
      layout = fit[["linear_design"]]
      expect_false(anyDuplicated(layout[["columns"]][["id"]]) > 0L)
      close(fit)
    }
  }

  # One-level nuisance factors pass through the same actual matrix producer.
  # Independent explicit group/trend columns define the fitted-value oracle.
  for (fml in list(y ~ x | g, y ~ x | g + h, y ~ x | g[x])) {
    fit = dbreg(fml, data = d, weights = "w", strategy = "compress")
    X = if (identical(deparse(fml), "y ~ x | g")) cbind(1, d[["x"]]) else
      if (identical(deparse(fml), "y ~ x | g + h")) model.matrix(~ x + h, d) else cbind(1, d[["x"]])
    oracle = lm.wfit(X, d[["y"]], d[["w"]])
    expect_equal(predict(fit, d), as.numeric(oracle[["fitted.values"]]), tolerance = 1e-8)
    expect_equal(predict(fit), predict(fit, fit[["data"]]), tolerance = 1e-8)
    expect_equal(df.residual(fit), n - if (identical(deparse(fml), "y ~ x | g + h")) 4 else 2)
    columns = fit[["linear_design"]][["columns"]]
    nuisance = which(columns[["role"]] == "fe_intercept")
    expect_true(any(vapply(columns[["recipe"]][nuisance], function(recipe) {
      identical(recipe[["kind"]], "constant") && identical(recipe[["variable"]], "g")
    }, logical(1))))
    close(fit)
  }

  # Retained R/SQL routes preserve missing placement, required sources and
  # caller resources. Generated-looking extras never replace raw categories.
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  d[["id"]] = seq_len(n)
  DBI::dbWriteTable(conn, "constant_source", d)
  for (retention in c("none", "prediction", "inference")) {
    fit = dbreg(y ~ 0 + cat + x, conn = conn, table = "constant_source",
      strategy = "compress", retain = retention)
    reference = lm.wfit(cbind(1, d[["x"]]), d[["y"]], rep(1, n))[["fitted.values"]]
    expect_equal(predict(fit, d), as.numeric(reference), tolerance = 1e-8)
    out = predict(fit, conn = conn, table = "constant_source", row_id = "id")
    actual = as.data.frame(out)
    expect_equal(actual[[out[["prediction_column"]]]][match(d[["id"]], actual[["id"]])],
      as.numeric(reference), tolerance = 1e-8)
    close(out)
    nd = d[1:4, , drop = FALSE]
    nd[["cat"]][2] = NA_character_
    expect_true(is.na(predict(fit, nd)[[2]]))
    nd[["cat"]] = NULL
    nd[["dbreg_linear_column1__"]] = 1
    expect_error(predict(fit, nd), "cat")
    close(fit)
    expect_true(DBI::dbIsValid(conn))
    expect_true(DBI::dbExistsTable(conn, "constant_source"))
  }

  # Producer controls: unused declared levels are not collapsed; numeric7
  # stays7; unaffected two-level matrix attributes and values are unchanged.
  control = data.frame(cat = factor(rep("one", 4), levels = c("one", "unused")), x = rep(7, 4))
  expect_equal(dbreg:::fitted_model_matrix(~ cat + x, control),
    model.matrix(~ cat + x, control, contrasts.arg = list(cat = contr.treatment(levels(control[["cat"]])))))
  expect_equal(dbreg:::fitted_model_matrix(~ 0 + x, control), model.matrix(~ 0 + x, control))
  control[["cat"]] = c(FALSE, TRUE, FALSE, TRUE)
  expect_equal(dbreg:::fitted_model_matrix(~ cat + x, control),
    model.matrix(~ cat + x, control, contrasts.arg = list(cat = contr.treatment(c("FALSE", "TRUE")))))

  # Structured constant leaves retain source identity through quoted names,
  # products and supported powers, without changing the public SQL compiler.
  names(d)[names(d) == "cat"] = "single category"
  for (fml in list(~ 0 + `single category`, ~ 0 + `single category`:h,
    ~ `single category`:I(x^2))) {
    X = dbreg:::fitted_model_matrix(fml, d)
    numeric_data = d
    numeric_data[["single category"]] = 1
    reference = model.matrix(fml, numeric_data)
    expect_equal(as.numeric(X), as.numeric(reference))
    expect_equal(attr(X, "assign"), attr(reference, "assign"))
    expect_identical(attr(X, "dbreg_constant_sources"), "single category")
    recipes = dbreg:::model_matrix_recipes(X, fml, d)
    expect_equal(vapply(recipes, dbreg:::linear_recipe_label, character(1)), colnames(X))
  }
})
