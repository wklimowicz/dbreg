# This independent oracle uses source matrices, never package names/term expansion.
local({
  set.seed(6231)
  n = 360L
  d = data.frame(x = rep(c(-1, 0, 1), 120L),
    f = rep(c("a", "b"), each = 180L), period = rep(1:3, each = 120L),
    z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n),
    w = rep(c(0.4, 1.2, 2), 120L), cluster = rep(1:18, each = 20L))
  d[["d1"]] = 0.8 * d[["z1"]] + 0.3 * d[["z3"]] + rnorm(n)
  d[["d2"]] = 0.7 * d[["z2"]] - 0.2 * d[["z3"]] + rnorm(n)
  X = cbind(1, d[["x"]], as.numeric(d[["f"]] == "b"),
    as.numeric(d[["period"]] == 2), as.numeric(d[["period"]] == 3),
    d[["d1"]], d[["d2"]])
  Z = cbind(X[, 1:5, drop = FALSE], d[["z1"]], d[["z2"]], d[["z3"]])
  d[["y"]] = drop(X %*% c(0.3, -0.4, 0.5, 0.2, -0.1, 0.7, -0.6)) + rnorm(n)
  w = d[["w"]]
  zwz = crossprod(Z, w * Z)
  xwz = crossprod(X, w * Z)
  P = xwz %*% solve(zwz)
  Ainv = solve(P %*% t(xwz))
  beta = Ainv %*% P %*% crossprod(Z, w * d[["y"]])
  residual = d[["y"]] - drop(X %*% beta)
  scores = (w * residual) * Z
  V_iid = Ainv * sum(w * residual^2) / (n - ncol(X))
  V_hc1 = Ainv %*% P %*% crossprod(scores) %*% t(P) %*% Ainv * n / (n - ncol(X))
  grouped = rowsum(scores, d[["cluster"]], reorder = FALSE)
  G = nrow(grouped)
  V_cluster = Ainv %*% P %*% crossprod(grouped) %*% t(P) %*% Ainv *
    G / (G - 1) * (n - 1) / (n - ncol(X))
  expected = list(iid = V_iid, hc1 = V_hc1, cluster = V_cluster)
  for (strategy in c("moments", "compress")) for (kind in names(expected)) {
    fit = dbreg::dbivreg(y ~ x + f + i(period, ref = 1) | d1 + d2 ~ z1 + z2 + z3,
      data = d, weights = "w", strategy = strategy,
      vcov = if (kind == "cluster") ~cluster else kind)
    expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
    expect_equal(stats::coef(fit, collin = TRUE), stats::coef(fit))
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(expected[[kind]]), tolerance = 1e-7)
    expect_equal(stats::predict(fit, d), drop(X %*% beta), tolerance = 1e-8)
    expect_equal(stats::nobs(fit), n)
    expect_equal(stats::df.residual(fit), n - ncol(X))
    expect_true(is.null(fit[["fitted_state"]]))
    expect_true(is.null(fit[["resource_state"]]))
    expect_true(is.null(fit[["data"]]))
    layout = fit[["linear_design"]]
    expect_equal(length(layout[["fitted_ids"]]), ncol(X))
    expect_equal(anyDuplicated(layout[["columns"]][["id"]]), 0L)
    expect_true(all(startsWith(layout[["columns"]][["id"]], "iv:")))
    expect_true(!any(startsWith(layout[["columns"]][["id"]], "iv:excluded:")))
    expect_equal(attr(stats::vcov(fit), "dbreg_coefficient_ids"), layout[["fitted_ids"]])
    # First-stage diagnostics use the full Z axis, not the structural X axis.
    for (j in 1:2) {
      outcome = d[[c("d1", "d2")[[j]]]]
      first_beta = solve(zwz, crossprod(Z, w * outcome))
      first_residual = outcome - drop(Z %*% first_beta)
      restricted_X = X[, 1:5, drop = FALSE]
      restricted_beta = solve(crossprod(restricted_X, w * restricted_X), crossprod(restricted_X, w * outcome))
      restricted_residual = outcome - drop(restricted_X %*% restricted_beta)
      rss = sum(w * first_residual^2)
      first_F = (sum(w * restricted_residual^2) - rss) / 3 / (rss / (n - ncol(Z)))
      expect_equal(fit[["diagnostics"]][["first_stage_f"]][[j]][["stat"]], first_F, tolerance = 1e-7)
      inverse = solve(zwz)
      first_scores = (w * first_residual) * Z
      first_V = if (kind == "iid") inverse * rss / (n - ncol(Z)) else if (kind == "hc1")
        inverse %*% crossprod(first_scores) %*% inverse * n / (n - ncol(Z)) else {
          sums = rowsum(first_scores, d[["cluster"]], reorder = FALSE)
          inverse %*% crossprod(sums) %*% inverse * G / (G - 1) * (n - 1) / (n - ncol(Z))
        }
      selected = 6:8
      first_wald = drop(t(first_beta[selected, , drop = FALSE]) %*%
        solve(first_V[selected, selected, drop = FALSE], first_beta[selected, , drop = FALSE])) / 3
      expect_equal(fit[["diagnostics"]][["first_stage_wald"]][[j]][["stat"]], first_wald, tolerance = 1e-7)
    }
    # Instruments are unnecessary in newdata once the structural equation is fit.
    raw = d[seq(1L, n, by = 17L), setdiff(names(d), c("z1", "z2", "z3")), drop = FALSE]
    expect_equal(stats::predict(fit, raw), drop(X[seq(1L, n, by = 17L), , drop = FALSE] %*% beta), tolerance = 1e-8)
  }
})

local({
  set.seed(6232)
  n = 320L
  d = data.frame(f = rep(c("a", "b"), 160L), fb = rnorm(n),
    p = rep(1:2, each = 160L), q = sample(rep(1:4, 80L)),
    z = rnorm(n), w = rep(c(0.5, 1.5), 160L))
  d[["endog"]] = d[["z"]] + d[["q"]] / 3 + rnorm(n)
  # Raw fb and the ordinary f='b' basis share a display label.
  X = cbind(1, d[["fb"]], as.numeric(d[["f"]] == "b"),
    as.numeric(d[["p"]] == 2), d[["endog"]])
  Z = cbind(X[, 1:4, drop = FALSE], d[["z"]],
    as.numeric(d[["q"]] == 2), as.numeric(d[["q"]] == 3), as.numeric(d[["q"]] == 4))
  d[["y"]] = drop(X %*% c(0.1, -0.5, 0.3, 0.2, 0.7)) + rnorm(n)
  w = d[["w"]]
  P = crossprod(X, w * Z) %*% solve(crossprod(Z, w * Z))
  Ainv = solve(P %*% crossprod(Z, w * X))
  beta = Ainv %*% P %*% crossprod(Z, w * d[["y"]])
  residual = d[["y"]] - drop(X %*% beta)
  V = Ainv %*% P %*% crossprod((w * residual) * Z) %*% t(P) %*% Ainv * n / (n - ncol(X))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbivreg(y ~ fb + f + i(p, ref = 1) | endog ~ z + i(q, ref = 1),
      data = d, strategy = strategy, weights = "w", vcov = "hc1")
    expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
    expect_equal(stats::coef(fit, collin = TRUE), stats::coef(fit))
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V), tolerance = 1e-7)
    rows = rev(seq_len(n))
    raw = d[rows, setdiff(names(d), c("q", "z")), drop = FALSE]
    raw[["dbreg_i1"]] = 917
    raw[["dbreg_internal_x_0001__"]] = -811
    expect_equal(stats::predict(fit, raw), drop(X[rows, , drop = FALSE] %*% beta), tolerance = 1e-8)
    expect_error(stats::predict(fit), "newdata is required")
    expect_error(dbreg::dbivreg(y ~ fb | endog ~ z, data = d, retain = "all"), "Unused arguments")
    expect_error(stats::vcov(fit, vcov = "hc1"), "requires fitting")
    expect_true(is.null(fit[["data"]]))
    expect_true(is.null(fit[["fitted_state"]]))
    layout = fit[["linear_design"]]
    expect_true(anyDuplicated(names(stats::coef(fit))) > 0L)
    expect_equal(anyDuplicated(layout[["columns"]][["id"]]), 0L)
    expect_true(all(c("iv:exog:", "iv:endog:") %in%
      unique(sub("t[0-9]+b[0-9]+$", "", layout[["columns"]][["id"]]))))
    broken = fit
    broken[["linear_design"]][["columns"]][["recipe"]][2L] = list(NULL)
    expect_error(stats::predict(broken, raw), "complete raw-data prediction adapter")
  }
})

# Numerical no-intercept IV remains a distinct emitted axis: no synthetic constant.
local({
  set.seed(6233)
  d = data.frame(x = rnorm(180L), z = rnorm(180L), z2 = rnorm(180L))
  d[["endog"]] = 0.7 * d[["z"]] + 0.2 * d[["z2"]] + rnorm(180L)
  d[["y"]] = -0.3 * d[["x"]] + 0.8 * d[["endog"]] + rnorm(180L)
  X = cbind(d[["x"]], d[["endog"]])
  Z = cbind(d[["x"]], d[["z"]], d[["z2"]])
  P = crossprod(X, Z) %*% solve(crossprod(Z))
  beta = solve(P %*% crossprod(Z, X), P %*% crossprod(Z, d[["y"]]))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbivreg(y ~ 0 + x | endog ~ z + z2, data = d, strategy = strategy)
    expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
    expect_equal(stats::coef(fit, collin = TRUE), stats::coef(fit))
    expect_equal(stats::predict(fit, d), drop(X %*% beta), tolerance = 1e-8)
    expect_true(!"iv:intercept" %in% fit[["linear_design"]][["fitted_ids"]])
    expect_error(stats::predict(fit), "newdata is required")
  }
})

# FE IV installs only structural inference axes; it does not activate learned FE prediction.
local({
  set.seed(6234)
  n = 240L
  d = data.frame(x = rnorm(n), z = rnorm(n), z2 = rnorm(n),
    g = rep(1:12, each = 20L), w = rep(c(0.5, 1, 1.5), 80L))
  d[["endog"]] = d[["z"]] + 0.3 * d[["z2"]] + rnorm(n)
  d[["y"]] = -0.3 * d[["x"]] + 0.8 * d[["endog"]] + d[["g"]] / 7 + rnorm(n)
  D = outer(d[["g"]], 1:12, `==`) * 1.0
  w = d[["w"]]
  raw_X = cbind(d[["x"]], d[["endog"]])
  raw_Z = cbind(d[["x"]], d[["z"]], d[["z2"]])
  within = function(value) value - D %*% solve(crossprod(D, w * D), crossprod(D, w * value))
  X = within(raw_X)
  Z = within(raw_Z)
  y = within(matrix(d[["y"]], ncol = 1L))
  P = crossprod(X, w * Z) %*% solve(crossprod(Z, w * Z))
  inverse = solve(P %*% crossprod(Z, w * X))
  beta = inverse %*% P %*% crossprod(Z, w * y)
  residual = y - X %*% beta
  V = inverse * sum(w * residual^2) / (n - ncol(X) - ncol(D))
  fit = dbreg::dbivreg(y ~ x | g | endog ~ z + z2, data = d, weights = "w", strategy = "demean")
  expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V), tolerance = 1e-7)
  expect_equal(stats::coef(fit, fe = TRUE), stats::coef(fit))
  expect_equal(stats::df.residual(fit), n - ncol(X) - ncol(D))
  expect_true(!dbreg:::native_linear_runtime(fit))
  expect_true(is.null(fit[["fitted_state"]]))
  expect_true(is.null(fit[["data"]]))
})

# Empty included-exogenous block preserves a genuine zero-row compiler result.
local({
  set.seed(6235)
  d = data.frame(z = rnorm(180L), z2 = rnorm(180L))
  d[["endog"]] = 0.8 * d[["z"]] - 0.2 * d[["z2"]] + rnorm(180L)
  d[["y"]] = 0.4 + 0.7 * d[["endog"]] + rnorm(180L)
  for (constant in c(TRUE, FALSE)) {
    X = if (constant) cbind(1, d[["endog"]]) else matrix(d[["endog"]], ncol = 1L)
    Z = if (constant) cbind(1, d[["z"]], d[["z2"]]) else cbind(d[["z"]], d[["z2"]])
    P = crossprod(X, Z) %*% solve(crossprod(Z))
    inverse = solve(P %*% crossprod(Z, X))
    beta = inverse %*% P %*% crossprod(Z, d[["y"]])
    residual = d[["y"]] - drop(X %*% beta)
    V = inverse * sum(residual^2) / (nrow(X) - ncol(X))
    fml = if (constant) y ~ 1 | endog ~ z + z2 else y ~ 0 | endog ~ z + z2
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbivreg(fml, data = d, strategy = strategy)
      expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V), tolerance = 1e-7)
      expect_equal(stats::predict(fit, d), drop(X %*% beta), tolerance = 1e-8)
      expect_equal("iv:intercept" %in% fit[["linear_design"]][["fitted_ids"]], constant)
      expect_true(!any(startsWith(fit[["linear_design"]][["columns"]][["id"]], "iv:exog:")))
    }
  }
  expect_error(dbreg::dbivreg(y ~ 1 | 0 ~ z, data = d), "endogenous regressor")
  expect_error(dbreg::dbivreg(y ~ 1 | endog ~ 0, data = d), "excluded instrument")
})

# A repeated indicator declaration in X and excluded Z remains the existing
# included-instrument overlap warning/removal, not a second saved dictionary.
local({
  set.seed(6236)
  d = data.frame(p = rep(1:3, 80L), z = rnorm(240L), z2 = rnorm(240L))
  d[["endog"]] = d[["z"]] + d[["z2"]] / 3 + rnorm(240L)
  X = cbind(1, as.numeric(d[["p"]] == 2), as.numeric(d[["p"]] == 3), d[["endog"]])
  Z = cbind(X[, 1:3, drop = FALSE], d[["z"]], d[["z2"]])
  d[["y"]] = drop(X %*% c(0.2, -0.3, 0.4, 0.7)) + rnorm(240L)
  P = crossprod(X, Z) %*% solve(crossprod(Z))
  beta = solve(P %*% crossprod(Z, X), P %*% crossprod(Z, d[["y"]]))
  for (strategy in c("moments", "compress")) {
    fit = NULL
    expect_warning({fit = dbreg::dbivreg(
      y ~ i(p, ref = 1) | endog ~ i(p, ref = 1) + z + z2, data = d, strategy = strategy)},
      "already valid instruments")
    expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
    expect_equal(stats::predict(fit, d), drop(X %*% beta), tolerance = 1e-8)
    expect_equal(length(fit[["i_specs"]]), 1L)
  }
  expect_error(dbreg::dbivreg(y ~ i(p, ref = 1) | i(p, ref = 1) ~ z + z2, data = d),
    "cannot also appear as included exogenous")
})

# Adapter metadata must not change the existing rank/identification refusals.
local({
  set.seed(6238)
  d = data.frame(z = rnorm(120L), z2 = rnorm(120L))
  d[["endog"]] = d[["z"]] + rnorm(120L)
  d[["endog2"]] = d[["z2"]] + rnorm(120L)
  d[["zcopy"]] = d[["z"]]
  d[["y"]] = d[["endog"]] + rnorm(120L)
  for (strategy in c("moments", "compress")) {
    expect_error(dbreg::dbivreg(y ~ 1 | endog + endog2 ~ z, data = d, strategy = strategy),
      "at least as many excluded instrument")
    expect_error(dbreg::dbivreg(y ~ 1 | endog ~ z + zcopy, data = d, strategy = strategy),
      "Instrument matrix is rank deficient")
  }
})

# Typed indicator dictionaries retain adjacent doubles in structural R recipes.
local({
  set.seed(6239)
  keys = 1 + c(0, 1, 2) * .Machine[["double.eps"]]
  d = data.frame(key = rep(keys, 80L), z = rnorm(240L), z2 = rnorm(240L))
  d[["endog"]] = d[["z"]] + 0.3 * d[["z2"]] + rnorm(240L)
  X = cbind(1, as.numeric(d[["key"]] == keys[[2L]]),
    as.numeric(d[["key"]] == keys[[3L]]), d[["endog"]])
  Z = cbind(X[, 1:3, drop = FALSE], d[["z"]], d[["z2"]])
  d[["y"]] = drop(X %*% c(0.2, -0.4, 0.6, 0.7)) + rnorm(240L)
  P = crossprod(X, Z) %*% solve(crossprod(Z))
  beta = solve(P %*% crossprod(Z, X), P %*% crossprod(Z, d[["y"]]))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbivreg(y ~ i(key, ref = 1) | endog ~ z + z2, data = d, strategy = strategy)
    expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
    expect_equal(stats::predict(fit, d), drop(X %*% beta), tolerance = 1e-8)
    unseen = d[1L, , drop = FALSE]
    unseen[["key"]] = 1 + 3 * .Machine[["double.eps"]]
    expect_error(stats::predict(fit, unseen), "unseen i.*levels")
  }
})
