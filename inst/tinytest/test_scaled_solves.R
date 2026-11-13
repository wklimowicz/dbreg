# Units must not determine numerical rank. B1a's original-scale AP floor
# remains necessary: normalization must not promote projection roundoff.
local({
  saved = options(dbreg.ap_tol = 1e-10, dbreg.ap_max_iter = 200L)
  on.exit(options(saved), add = TRUE)
  set.seed(711)
  n = 240L
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(1:9, n, TRUE),
    h = sample(1:7, n, TRUE), w = runif(n, .4, 2), cl = rep(1:12, each = 20))
  d[["endog"]] = .8*d[["z"]] + .2*d[["x"]] + rnorm(n)
  d[["y"]] = .7 + .4*d[["x"]] + .8*d[["endog"]] + sin(d[["g"]]) + rnorm(n)
  for (iv in c(FALSE, TRUE)) {
    for (path in c("moments", "compress", "demean", "ap")) {
      strategy = if (path == "ap") "demean" else path
      formula = if (iv) {
        if (path == "ap") y ~ x | g + h | endog ~ z else if (path == "demean") y ~ x | g | endog ~ z else y ~ x | endog ~ z
      } else if (path == "ap") y ~ x | g + h else if (path == "demean") y ~ x | g else y ~ x
      for (kind in list("iid", "hc1", ~cl)) {
        ref = fixest::feols(formula, d, weights = ~w, vcov = kind, fixef.rm = "none",
          ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
        b = coef(ref); names(b) = sub("^fit_", "", names(b))
        V = vcov(ref); dimnames(V) = lapply(dimnames(V), function(x) sub("^fit_", "", x))
        # AP retains its absolute convergence gate; test solver scaling within
        # that gate, and wider units on direct and analytic-demeaning paths.
        for (unit in if (path == "ap") c(1e-4, 1e4) else c(1e-12, 1e12)) {
          scaled = d
          scaled[["x"]] = scaled[["x"]]*unit
          scaled[["endog"]] = scaled[["endog"]]/unit
          scaled[["z"]] = scaled[["z"]]*unit
          estimator = if (iv) dbivreg else dbreg
          fit = tryCatch(estimator(formula, data = scaled, weights = "w", vcov = kind, strategy = strategy), error = function(e) e)
          expect_false(inherits(fit, "error"), info = paste(iv, path, deparse(kind), unit, if (inherits(fit, "error")) conditionMessage(fit) else ""))
          if (inherits(fit, "error")) next
          terms = names(b)
          expect_true(setequal(names(coef(fit)), terms), info = paste(iv, path, deparse(kind), unit, paste(names(coef(fit)), collapse = ","), "oracle:", paste(terms, collapse = ",")))
          if (!setequal(names(coef(fit)), terms)) next
          conversion = stats::setNames(rep(1, length(terms)), terms)
          conversion["x"] = unit
          if (iv) conversion["endog"] = 1/unit
          beta = coef(fit)[terms]*conversion
          covariance = as.matrix(vcov(fit))[terms, terms, drop = FALSE]*outer(conversion, conversion)
          expect_true(max(abs(beta-b)) < 1e-7)
          expect_true(max(abs(covariance-V)) < 1e-7)
        }
      }
    }
  }
  # A dependent column stays dependent at extreme but finite scales.
  X = cbind("(Intercept)" = 1, x = d[["x"]]*1e-12, copy = 2*d[["x"]]*1e12)
  S = crossprod(cbind(X, y = d[["y"]]))
  fit = dbreg:::fit_ols(S, "y", colnames(X), drop_collinear = TRUE)
  expect_identical(fit[["dropped"]], "copy")
  expect_identical(fit[["x"]], c("(Intercept)", "x"))
})

# Equilibrate both vector and matrix right-hand sides, then restore units.
local({
  set.seed(19)
  X = cbind(1, rnorm(30), rnorm(30))
  A = crossprod(X)
  B = cbind(rnorm(3), rnorm(3))
  units = c(1, 1e-12, 1e12)
  scaled = A*outer(units, units)
  fit = dbreg:::solve_with_fallback(scaled, B*units)
  expect_equal(unname(fit[["betahat"]]*units), unname(solve(A, B)), tolerance = 1e-12)
  expect_equal(unname(fit[["XtX_inv"]]*outer(units, units)), unname(solve(A)), tolerance = 1e-12)
  vector_fit = dbreg:::solve_with_fallback(scaled, B[, 1L]*units)
  expect_equal(as.numeric(vector_fit[["betahat"]])*units, as.numeric(solve(A, B[, 1L])), tolerance = 1e-12)
  # The QR fallback still handles a nonsingular system without a Cholesky factor.
  indefinite = matrix(c(1, 2, 2, 1), 2)
  fallback = dbreg:::solve_with_fallback(indefinite, c(1, 3))
  expect_equal(as.numeric(fallback[["betahat"]]), as.numeric(solve(indefinite, c(1, 3))), tolerance = 1e-12)
  expect_equal(fallback[["XtX_inv"]], solve(indefinite), tolerance = 1e-12)
  with_zero = cbind(X, zero = 0)
  collin = dbreg:::detect_collinearity(crossprod(with_zero), crossprod(with_zero, rnorm(30)))
  expect_identical(collin[["drop_names"]], "zero")
})
