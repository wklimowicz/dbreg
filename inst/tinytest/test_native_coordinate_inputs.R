# Preserve raw-source validation and the explicit fitted-owned cell adapter.
local({
  set.seed(11012)
  d = data.frame(x = rep(1:3, 40), z = stats::rnorm(120))
  d[["y"]] = 0.3 * d[["x"]] + (d[["x"]] == 2) + 0.7 * d[["z"]] + stats::rnorm(120, sd = 0.1)
  for (formula in list(y ~ z + factor(x), y ~ x + z + factor(x))) {
    reference = stats::lm(formula, d)
    reference_beta = stats::coef(reference)
    reference_beta[is.na(reference_beta)] = 0
    rows = data.frame(x = c(1, 2, 3), z = c(0.2, -0.3, 0.4))
    reference_X = stats::model.matrix(stats::delete.response(stats::terms(reference)), rows,
      contrasts.arg = reference[["contrasts"]], xlev = reference[["xlevels"]])
    expected = as.numeric(reference_X %*% as.numeric(reference_beta))
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbreg(formula, data = d, strategy = strategy)
      expect_equal(stats::predict(fit, rows), expected, tolerance = 1e-9)
      bad = data.frame(x = 99, z = 0.3)
      bad[["factor(x)"]] = 1 # Generated-looking extras cannot override raw99.
      expect_error(stats::predict(fit, bad), "unseen levels.*factor")
      filtered = bad
      filtered[["z"]] = NA_real_
      expect_identical(stats::predict(fit, filtered), NA_real_)
      close(fit)
    }
  }
  one = d
  one[["x"]] = 2
  reference = stats::lm(y ~ z, one)
  rows = data.frame(x = c(2, 2), z = c(0.2, -0.3))
  expected = as.numeric(cbind(1, rows[["z"]]) %*% as.numeric(stats::coef(reference)))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbreg(y ~ factor(x) + z, data = one, strategy = strategy)
    expect_equal(stats::predict(fit, rows), expected, tolerance = 1e-9)
    expect_error(stats::predict(fit, data.frame(x = 99, z = 0)), "unseen levels.*factor")
    close(fit)
  }
})

local({
  # Every cell mean is exactly linear in an independently constructed matrix.
  # Distinct prepared cells correspond to unique rows of that matrix. Thus
  # sorted cell predictions have an independent raw-design oracle, without
  # using saved SQL aliases, fitted labels or package coefficients as truth.
  d = expand.grid(x = 1:4, f = c("a", "b", "c"), z = 1:3, replicate = 1:3,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  cases = list(
    list(formula = y ~ I(x^2), X = cbind(1, d[["x"]]^2)),
    list(formula = y ~ i(f, ref = "a"), X = cbind(1, as.numeric(d[["f"]] == "b"), as.numeric(d[["f"]] == "c"))),
    list(formula = y ~ factor(f):z, X = cbind(1, as.numeric(d[["f"]] == "a") * d[["z"]],
      as.numeric(d[["f"]] == "b") * d[["z"]], as.numeric(d[["f"]] == "c") * d[["z"]]))
  )
  for (case in cases) {
    X = case[["X"]]
    beta = seq_len(ncol(X)) + sqrt(seq_len(ncol(X))) / 10
    d[["y"]] = as.numeric(X %*% beta) + 0.1 * (d[["replicate"]] - 2)
    expected_cells = sort(as.numeric(as.matrix(unique(as.data.frame(X))) %*% beta))
    fit = dbreg::dbreg(case[["formula"]], data = d, strategy = "compress")
    cells = fit[["data"]]
    expect_equal(sort(as.numeric(stats::predict(fit))), expected_cells, tolerance = 1e-8)
    expect_equal(stats::predict(fit, cells), stats::predict(fit), tolerance = 1e-10)
    # A changed arbitrary frame containing prepared aliases is not fitted-owned.
    # It still must supply the raw formula sources when those are absent.
    if (!all(all.vars(stats::formula(fit[["fml"]], lhs = 0, rhs = 1)) %in% names(cells))) {
      arbitrary = cells
      arbitrary[["extra"]] = 1
      expect_error(stats::predict(fit, arbitrary), "missing required columns")
    }
    close(fit)
  }
})

local({
  # Preserve exact numeric comparison values before display factor coercion,
  # while caller-supplied factor labels retain their existing label semantics.
  set.seed(11014)
  values = c(0.12345678901234567, 0.76543210987654321, 1.2345678901234567)
  d = data.frame(x = stats::rnorm(120), category = rep(values, 40))
  X = cbind(1, d[["x"]], as.numeric(d[["category"]] == values[[2L]]), as.numeric(d[["category"]] == values[[3L]]))
  d[["y"]] = as.numeric(X %*% c(1, 0.5, -0.4, 0.3)) + stats::rnorm(120, sd = 0.1)
  beta = as.numeric(qr.solve(X, d[["y"]]))
  rows = data.frame(x = c(-0.2, 0.4, 0.6), category = values)
  expected = as.numeric(cbind(1, rows[["x"]], c(0, 1, 0), c(0, 0, 1)) %*% beta)
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbreg(y ~ x + category, data = d, strategy = strategy,
      factor_refs = list(category = values[[1L]]))
    expect_equal(stats::predict(fit, rows), expected, tolerance = 1e-10)
    labelled = rows
    labelled[["category"]] = factor(as.character(values), levels = as.character(values))
    expect_equal(stats::predict(fit, labelled), expected, tolerance = 1e-10)
    close(fit)
  }
})
