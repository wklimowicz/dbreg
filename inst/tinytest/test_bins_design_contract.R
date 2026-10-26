# Plot levels and their uncertainty cannot depend on automatic strategy choice.
local({
  set.seed(1203)
  d = data.frame(x = stats::runif(480), g = rep(1:160, each = 3))
  d[["y"]] = 10 + 3 * d[["x"]] + sin(d[["g"]]) + stats::rnorm(nrow(d))
  breaks = seq(0, 1, length.out = 5)
  d[["b"]] = cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE)
  indicators = stats::model.matrix(~b, d)
  for (degree in 0:1) {
    u = (d[["x"]] - as.numeric(tapply(d[["x"]], d[["b"]], min))[as.integer(d[["b"]])]) /
      (as.numeric(tapply(d[["x"]], d[["b"]], max)) - as.numeric(tapply(d[["x"]], d[["b"]], min)))[as.integer(d[["b"]])]
    slopes = if (degree) stats::model.matrix(~b - 1, d) * u else NULL
    X = cbind(indicators, slopes, stats::model.matrix(~factor(g), d)[, -1, drop = FALSE])
    fit = stats::lm.fit(X, d[["y"]])
    bread = solve(crossprod(X))
    Z = matrix(0, 4L, ncol(X))
    Z[, seq_len(4L)] = cbind(1, diag(4L)[, -1, drop = FALSE])
    if (degree) Z[, 4L + seq_len(4L)] = diag(as.numeric(tapply(u, d[["b"]], mean)))
    expected = as.numeric(Z %*% fit[["coefficients"]])
    for (kind in c("iid", "hc1")) {
      V = if (kind == "iid") bread * sum(fit[["residuals"]]^2) / (nrow(X) - ncol(X)) else {
        scores = X * fit[["residuals"]]
        nrow(X) / (nrow(X) - ncol(X)) * bread %*% crossprod(scores) %*% bread
      }
      standard_errors = sqrt(diag(Z %*% V %*% t(Z)))
      models = lapply(c("auto", "compress"), function(strategy) dbbinsreg(y ~ x | g, data = d,
        nbins = 4, binspos = breaks, points = c(degree, 0), strategy = strategy, vcov = kind, plot = FALSE))
      for (model in models) {
        expect_true(all(is.finite(model[["points"]][["fit"]])))
        expect_equal(model[["points"]][["fit"]], expected, tolerance = 1e-8)
        expect_equal(model[["points"]][["se"]], standard_errors, tolerance = 1e-8)
      }
      expect_equal(models[[1]][["points"]], models[[2]][["points"]], tolerance = 1e-10)
    }
  }
})

# Names used by the binning algorithm must not shadow outcomes/controls/keys.
local({
  set.seed(1205)
  d = data.frame(x = stats::runif(240), z = stats::rnorm(240), g = rep(1:12, each = 20))
  d[["y"]] = 1 + 3 * d[["x"]] + 2 * d[["z"]] + stats::rnorm(nrow(d))
  for (degree in 0:1) {
    baseline = dbbinsreg(y ~ x + z, data = d, nbins = 4, points = c(degree, 0), plot = FALSE)
    for (name in c("bin", "bin1", "bin2", "u_1", "u2_2", "x1_spline")) {
      renamed = d
      names(renamed)[names(renamed) == "z"] = name
      actual = dbbinsreg(stats::reformulate(c("x", name), response = "y"), data = renamed,
        nbins = 4, points = c(degree, 0), plot = FALSE)
      expect_equal(actual[["points"]], baseline[["points"]], tolerance = 1e-8)
      expect_identical(anyDuplicated(rownames(actual[["model"]][["coeftable"]])), 0L)
    }
    # The running variable and outcome can themselves use reserved names.
    renamed = d
    names(renamed)[names(renamed) == "x"] = "bin"
    names(renamed)[names(renamed) == "y"] = "bin2"
    actual = dbbinsreg(bin2 ~ bin + z, data = renamed, nbins = 4, points = c(degree, 0), plot = FALSE)
    expect_equal(actual[["points"]], baseline[["points"]], tolerance = 1e-8)
  }
  # Returned compressed models keep a usable design after source restoration.
  for (name in c("bin2", "h")) {
    renamed = d
    names(renamed)[names(renamed) == "z"] = name
    model = dbbinsreg(stats::reformulate(c("x", name), response = "y"), data = renamed,
      nbins = 4, strategy = "compress", plot = FALSE)[["model"]]
    values = predict(model)
    data = model[["data"]]
    formula = stats::formula(model[["fml"]])
    formula[[2L]] = as.name("mean_Y")
    X = stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
    reference = stats::lm.wfit(X, data[["mean_Y"]], data[["n"]])
    expect_equal(as.numeric(values), as.numeric(X %*% reference[["coefficients"]]), tolerance = 1e-8)
  }
  # Exercise the constrained spline path as well as polynomial columns.
  baseline = dbbinsreg(y ~ x + z, data = d, nbins = 4, points = c(1, 1), plot = FALSE)
  renamed = d
  names(renamed)[names(renamed) == "z"] = "x1_spline"
  actual = dbbinsreg(y ~ x + x1_spline, data = renamed, nbins = 4, points = c(1, 1), plot = FALSE)
  expect_equal(actual[["points"]], baseline[["points"]], tolerance = 1e-8)
  expect_identical(anyDuplicated(rownames(actual[["model"]][["coeftable"]])), 0L)
})
