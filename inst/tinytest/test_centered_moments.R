# Translation and a common weight multiplier must preserve identified slopes,
# inference and residual variation. Oracles work in ordinary coordinates.
local({
  saved = options(dbreg.ap_max_iter = 200L)
  on.exit(options(saved), add = TRUE)
  set.seed(193)
  n = 180L
  d = data.frame(x = rnorm(n), z = rnorm(n), z2 = rnorm(n),
    w = runif(n, .4, 2), g = rep(1:9, each = 20), cl = rep(1:15, 12))
  d[["endog"]] = .7*d[["z"]]+.4*d[["z2"]]+.2*d[["x"]]+rnorm(n)
  d[["y"]] = .7+.4*d[["x"]]+.9*d[["endog"]]+sin(d[["g"]])+rnorm(n)
  shifts = c(x = 1e8, z = 3e8, z2 = -1e8, endog = -2e8, y = 1e8)
  for (iv in c(FALSE, TRUE)) for (strategy in c("moments", "compress", "demean")) {
    formula = if (iv) {
      if (strategy == "demean") y ~ x | g | endog ~ z + z2 else y ~ x | endog ~ z + z2
    } else if (strategy == "demean") y ~ x + z | g else y ~ x + z
    # Demeaned FE frames use the stored weights; extreme weight scales are out of scope there.
    for (scale in if (strategy == "demean") list(NULL, 1) else list(NULL, 1, 1e-200, 1e200)) for (kind in list("iid", "hc1", ~cl)) {
      ref = fixest::feols(formula, d, weights = if (!is.null(scale)) ~w else NULL,
        vcov = kind, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
      b = coef(ref); names(b) = sub("^fit_", "", names(b))
      V = vcov(ref); dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
      shifted = d
      for (name in names(shifts)) shifted[[name]] = shifted[[name]]+shifts[[name]]
      if (!is.null(scale)) shifted[["w"]] = shifted[["w"]]*scale
      estimator = if (iv) dbivreg else dbreg
      fit = tryCatch(estimator(formula, data = shifted, weights = if (!is.null(scale)) "w" else NULL,
        vcov = kind, strategy = strategy), error = function(e) e)
      info = paste(iv, strategy, if (is.null(scale)) "unweighted" else scale, deparse(kind))
      expect_false(inherits(fit, "error"), info = paste(info, if (inherits(fit, "error")) conditionMessage(fit) else ""))
      if (inherits(fit, "error")) next
      expect_true(setequal(names(coef(fit)), names(b)), info = info)
      if (!setequal(names(coef(fit)), names(b))) next
      terms = names(b)
      actual = coef(fit)[terms]
      A = diag(length(terms)); dimnames(A) = list(terms, terms)
      if ("(Intercept)" %in% terms) {
        slopes = setdiff(terms, "(Intercept)")
        A["(Intercept)", slopes] = -shifts[slopes]
        actual["(Intercept)"] = actual["(Intercept)"]-shifts[["y"]]+sum(shifts[slopes]*actual[slopes])
      }
      expect_true(max(abs(actual-b)) < 5e-6, info = info)
      expected = A %*% V %*% t(A)
      covariance = as.matrix(vcov(fit))[terms, terms, drop = FALSE]
      expect_true(all(abs(covariance-expected) < 1e-8+5e-6*abs(expected)), info = info)
      mass = if (is.null(scale)) 1 else scale
      expected_rss = sum(if (is.null(scale)) residuals(ref)^2 else d[["w"]]*residuals(ref)^2)
      expect_equal(as.numeric(attr(vcov(fit), "rss"))/mass, expected_rss, tolerance = 5e-6, info = info)
      expect_equal(nobs(fit), n, info = info)
      if (iv) {
        # Diagnostics must use moment coordinates too, including first stages.
        ordinary = estimator(formula, data = d, weights = if (!is.null(scale)) "w" else NULL,
          vcov = kind, strategy = strategy)
        expect_equal(fit[["diagnostics"]][["first_stage_f"]], ordinary[["diagnostics"]][["first_stage_f"]], tolerance = 5e-6, info = info)
        expect_equal(fit[["diagnostics"]][["first_stage_wald"]], ordinary[["diagnostics"]][["first_stage_wald"]], tolerance = 5e-6, info = info)
        expect_true(is.finite(fit[["diagnostics"]][["overid"]][["stat"]]), info = info)
        expect_equal(fit[["diagnostics"]][["overid"]], ordinary[["diagnostics"]][["overid"]], tolerance = 5e-6, info = info)
      }
    }
  }
})

local({
  set.seed(209)
  d = data.frame(x = rnorm(120), z = rnorm(120), w = runif(120, .4, 2),
    g = rep(1:6, each = 20), cl = rep(1:12, 10), f = rep(c("a", "b"), 60))
  d[["y"]] = .7+.4*d[["x"]]+.2*d[["z"]]+sin(d[["g"]])+rnorm(120)
  # No-intercept fits still normalize weights without changing their model.
  # Explicit intercept-only compression has a single cell and no grouping keys.
  for (formula in list(y ~ 0+x+z, y ~ 1, y ~ x+f | g)) {
    for (strategy in c("moments", "compress")) {
      if (length(formula) && identical(formula, y ~ x+f | g) && strategy == "moments") next
      for (kind in list("iid", "hc1", ~cl)) for (scale in c(1e-200, 1e200)) {
        scaled = d; scaled[["w"]] = scaled[["w"]]*scale
        fit = tryCatch(dbreg(formula, data = scaled, weights = "w", vcov = kind, strategy = strategy), error = function(e) e)
        expect_false(inherits(fit, "error"), info = paste(deparse(formula), strategy, deparse(kind), scale, if(inherits(fit,"error"))conditionMessage(fit)else""))
        if (inherits(fit, "error")) next
        ref = fixest::feols(formula, d, weights = ~w, vcov = kind, fixef.rm = "none",
          ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
        expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-7)
        expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-7)
      }
    }
  }
  # Force the cluster-by-cell route as well as the SQL literal-score route.
  for (iv in c(FALSE, TRUE)) for (scale in c(1e-200, 1e200)) {
    formula = if (iv) y ~ x | z ~ g else y ~ x+z
    estimator = if (iv) dbivreg else dbreg
    shifted = d; shifted[["x"]] = shifted[["x"]]+1e8
    shifted[["y"]] = shifted[["y"]]+1e8; shifted[["w"]] = shifted[["w"]]*scale
    a = tryCatch(estimator(formula, data = shifted, weights = "w", vcov = ~cl, strategy = "compress"), error = function(e) e)
    b = tryCatch(estimator(formula, data = shifted, weights = "w", vcov = ~cl, strategy = "compress", drop_missings = FALSE), error = function(e) e)
    expect_false(inherits(a, "error"))
    expect_false(inherits(b, "error"))
    if (inherits(a, "error") || inherits(b, "error")) next
    expect_equal(as.numeric(coef(a)), as.numeric(coef(b)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(a)), as.numeric(vcov(b)), tolerance = 1e-7)
  }
  # Internal location/scale fields must not collide with actual predictors.
  d[["dbreg_cell_center__"]] = d[["x"]]+1e8
  d[["dbreg_moment_center_1"]] = d[["z"]]-1e8
  d[["y"]] = d[["y"]]+1e8
  for (strategy in c("moments", "compress")) {
    fit = dbreg(y ~ dbreg_cell_center__+dbreg_moment_center_1, data = d, strategy = strategy, weights = "w", vcov = "hc1")
    ref = dbreg(y ~ x+z, data = transform(d, y = y-1e8), strategy = strategy, weights = "w", vcov = "hc1")
    expect_equal(as.numeric(coef(fit)[-1]), as.numeric(coef(ref)[-1]), tolerance = 1e-6)
    expect_equal(as.matrix(vcov(fit))[-1,-1], unname(as.matrix(vcov(ref))[-1,-1]), check.attributes = FALSE, tolerance = 1e-6)
  }
})

local({
  d = data.frame(x = rep(1:3, each = 20), y = 1e8+sin(1:60), w = 1+(1:60)%%3, cl = rep(1:10, 6))
  ref = fixest::feols(y ~ 1, d, weights = ~w, vcov = ~cl)
  for (ratio in c(.01, 1)) {
    fit = tryCatch(dbreg(y ~ 1, data = d, weights = "w", vcov = ~cl, compress_ratio = ratio), error = function(e) e)
    expect_false(inherits(fit, "error"))
    if (inherits(fit, "error")) next
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-12)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-6)
    expect_equal(nobs(fit), 60)
  }
  data = dbreg(y ~ x, data = d, weights = "w", strategy = "compress", data_only = TRUE)
  data = data[order(data[["x"]]), , drop = FALSE]
  expect_false(any(grepl("dbreg_cell_(center|scale)", names(data))))
  expected_w = as.numeric(rowsum(d[["w"]], d[["x"]]))
  expected_wy = as.numeric(rowsum(d[["w"]]*d[["y"]], d[["x"]]))
  expect_equal(data[["sum_w"]], expected_w)
  expect_equal(data[["sum_wy"]], expected_wy, tolerance = 1e-12)
  expect_equal(data[["mean_Y"]], expected_wy/expected_w, tolerance = 1e-12)
  expect_error(dbreg(y ~ 0, data = d), "No regressors")
})

# Nonidentity score maps require a transpose; raw weight mass must cancel.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  d = data.frame(x = (1:10)/10, e = sin(1:10), w = .5+(1:10)/10)
  X = cbind("(Intercept)" = 1, x = d[["x"]])
  expected = crossprod(X*(d[["e"]]*d[["w"]]))
  inputs = list2env(list(conn = conn, vcov_type_req = "hc1", cluster_var = NULL))
  frame = list(cte_sql = "WITH base AS (SELECT * FROM score_scale)", cte_name = "base", weights_expr = "w")
  for (scale in c(1e-200, 1e200)) {
    scaled = d; scaled[["w"]] = scaled[["w"]]*scale
    DBI::dbWriteTable(conn, "score_scale", scaled, overwrite = TRUE)
    frame[["weight_scale"]] = scale
    meat = dbreg:::compute_score_meat(inputs, frame, c(model = "e"), "x", score_intercept = TRUE)[[1L]]
    expect_identical(attr(meat, "weight_scale"), scale)
    expect_equal(as.numeric(meat), as.numeric(expected), tolerance = 1e-12)
    B = matrix(c(1, .4, 0, 2), 2)
    variance = dbreg:::compute_vcov(vcov_type = "hc1", XtX_inv = diag(2)/scale,
      sandwich_bread = B/scale, rss = 10*scale, df_res = 8, nobs_orig = 10, meat = meat)
    expect_equal(as.numeric(variance), as.numeric(10/8*B%*%expected%*%t(B)), tolerance = 1e-12)
  }
})

# Stabilizing moments must not leave compressed-cluster matching lossy.
local({
  set.seed(238)
  d = data.frame(x = 1e9+1e-5*rnorm(120), w = runif(120, .4, 2), cl = rep(1:12, 10))
  d[["u"]] = (d[["x"]]-1e9)/1e-5
  d[["y"]] = .4*d[["u"]]+rnorm(120)
  expect_true(length(unique(d[["x"]])) > length(unique(as.character(d[["x"]]))))
  ref = fixest::feols(y ~ u, d, weights = ~w, vcov = ~cl)
  for (drop in c(TRUE, FALSE)) {
    fit = dbreg(y ~ x, data = d, weights = "w", vcov = ~cl, strategy = "compress", drop_missings = drop)
    expect_equal(as.numeric(coef(fit)[["x"]])*1e-5, as.numeric(coef(ref)[["u"]]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)["x", "x"])*1e-10, as.numeric(vcov(ref)["u", "u"]), tolerance = 1e-7)
  }
  for (keys in list(data.frame(f = c("a\rb", "a"), h = c("c", "b\rc")),
                    data.frame(f = c(NA_character_, "NA"), h = c("c", "c")))) {
    cells = keys[c(1, 2, 1, 2), , drop = FALSE]
    cells[["cl"]] = c(1, 1, 2, 2)
    cells[["sum_w_gc"]] = 1
    cells[["sum_wy_gc"]] = 1:4
    X = diag(2); colnames(X) = c("A", "B")
    meat = dbreg:::cluster_meat_from_cells(cells, keys, c("f", "h"), "cl", X, c(.2, .4))
    expect_equal(as.numeric(meat), as.numeric(crossprod(matrix(c(.8, 2.8, 1.6, 3.6), 2))), tolerance = 1e-12)
    expect_identical(attr(meat, "n_clusters"), 2L)
  }
})

# Numeric FE labels become factors for the cell design, and must still match
# the same physical keys returned by the cluster query (e.g. double 0.1).
local({
  set.seed(487)
  d = data.frame(x = rnorm(120), g = rep((1:6)/10, each = 20), cl = rep(1:12, 10), w = runif(120, .4, 2))
  d[["y"]] = .7+.4*d[["x"]]+sin(d[["g"]])+rnorm(120)
  ref = fixest::feols(y ~ x | g, d, weights = ~w, vcov = ~cl, fixef.rm = "none",
    ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
  for (drop in c(TRUE, FALSE)) {
    fit = dbreg(y ~ x | g, data = d, weights = "w", vcov = ~cl, strategy = "compress", drop_missings = drop)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-7)
  }
})

# Sparse and dense design centering must match dense sweep() and preserve dimnames.
local({
  d = data.frame(x1 = c(1, 2, 3, 4), x2 = c(10, 20, 30, 40), g = factor(c("a", "b", "a", "b")))
  X_sp = Matrix::sparse.model.matrix(~ x1 + x2 + g, d)
  X_dn = as.matrix(X_sp)
  for (center in list(
    c(0, 2.5, 25, 0),
    c(0, 0, 0, 0),
    c(1, 2, 3, 4)
  )) {
    C_sp = dbreg:::center_compressed_design(X_sp, center)
    C_dn = dbreg:::center_compressed_design(X_dn, center)
    ref = sweep(X_dn, 2L, center, "-")
    expect_true(inherits(C_sp, "Matrix"))
    expect_equal(as.matrix(C_sp), ref, tolerance = 1e-15)
    expect_equal(C_dn, ref, tolerance = 1e-15)
    expect_identical(dimnames(C_sp), dimnames(X_sp))
  }
  empty_sp = X_sp[0, , drop = FALSE]
  expect_identical(dim(dbreg:::center_compressed_design(empty_sp, c(0, 1, 2, 0))), c(0L, ncol(X_sp)))
})
