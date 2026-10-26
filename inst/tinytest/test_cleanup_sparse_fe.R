local({
  graph = function(big, small, components = 1L, mobility = 0.001, scale = 1) {
    i = rep(seq_len(big), 2L)
    home = (seq_len(big) - 1L) %% small + 1L
    j = c(home, (home - 1L + components) %% small + 1L)
    h = c(rep(1, big), rep(mobility, length.out = big)) * scale
    info = list(list(id = seq_len(big), h = as.numeric(rowsum(h, i, reorder = TRUE))),
      list(id = seq_len(small), h = as.numeric(rowsum(h, j, reorder = TRUE))))
    pair = list(i = 1L, j = 2L, ir = i, jr = j, h = h)
    list(sys = dbreg:::fe_system(info, list(pair)), i = i, j = j, h = h)
  }
  set.seed(991)
  for (components in c(1L, 3L)) {
    fixture = graph(1000L, 600L, components)
    sys = fixture[["sys"]]
    cache = new.env(parent = emptyenv())
    solver = dbreg:::fe_solver(sys, cache, retain_null = TRUE)
    expect_identical(solver[["method"]], "sparse_cholesky")
    truth = matrix(rnorm(sys[["n"]] * 2L), sys[["n"]], 2L)
    rhs = dbreg:::fe_multiply(sys, truth)
    solution = solver[["solve"]](solver[["project"]](rhs))
    expect_equal(dbreg:::fe_multiply(sys, solution), rhs, tolerance = 1e-10)
    fitted = solution[fixture[["i"]], , drop = FALSE] + solution[1000L + fixture[["j"]], , drop = FALSE]
    expected = truth[fixture[["i"]], , drop = FALSE] + truth[1000L + fixture[["j"]], , drop = FALSE]
    expect_equal(fitted, expected, tolerance = 1e-8)
    null = solver[["null_basis"]]
    expect_equal(ncol(null), components)
    expect_equal(crossprod(null), diag(components), tolerance = 1e-12)
    expect_equal(dbreg:::fe_multiply(sys, null), matrix(0, sys[["n"]], components), tolerance = 1e-12)
    noise = matrix(rnorm(sys[["n"]] * 2L), sys[["n"]], 2L)
    projected = solver[["project"]](noise)
    expect_equal(crossprod(null, projected), matrix(0, components, 2L), tolerance = 1e-12)
    permutation = cache[["sparse_factor"]]@perm
    next_fixture = graph(1000L, 600L, components, mobility = 0.002)
    next_sys = next_fixture[["sys"]]
    next_solver = dbreg:::fe_solver(next_sys, cache)
    next_rhs = dbreg:::fe_multiply(next_sys, truth)
    expect_equal(dbreg:::fe_multiply(next_sys, next_solver[["solve"]](next_rhs)), next_rhs, tolerance = 1e-10)
    expect_identical(cache[["sparse_factor"]]@perm, permutation)
  }
  # Numerical zeros can split the information graph during IRLS even if
  # the CSC storage pattern is unchanged. Rebuild the component pins then.
  fixture = graph(1000L, 600L, mobility = 1)
  cache = new.env(parent = emptyenv())
  dbreg:::fe_solver(fixture[["sys"]], cache)
  home = (seq_len(1000L) - 1L) %% 600L + 1L
  mobility = as.numeric(!home %in% c(1L, 301L))
  split = graph(1000L, 600L, mobility = mobility)
  solver = dbreg:::fe_solver(split[["sys"]], cache)
  expect_equal(length(unique(cache[["sparse_components"]])), 2L)
  truth = matrix(rnorm(1600L), 1600L, 1L)
  rhs = dbreg:::fe_multiply(split[["sys"]], truth)
  expect_equal(dbreg:::fe_multiply(split[["sys"]], solver[["solve"]](rhs)), rhs, tolerance = 1e-10)
  # Divide stored nonzeros by the group weights directly: constructing their
  # reciprocals overflows for valid subnormal information weights.
  tiny = graph(1000L, 600L, mobility = 1, scale = 1e-310)
  solver = dbreg:::fe_solver(tiny[["sys"]], new.env(parent = emptyenv()))
  rhs = dbreg:::fe_multiply(tiny[["sys"]], truth)
  solution = solver[["solve"]](solver[["project"]](rhs))
  expect_true(all(is.finite(solution)))
  expect_equal(dbreg:::fe_multiply(tiny[["sys"]], solution) / 1e-310, rhs / 1e-310, tolerance = 1e-9)
  # The former dense cross block exceeded the default 512 MiB planning cap.
  fixture = graph(12500L, 1563L)
  sys = fixture[["sys"]]
  expect_true(inherits(sys[["Cm"]], "sparseMatrix"))
  expect_true(as.numeric(object.size(sys[["Cm"]])) < 3 * 1024^2)
  expect_true(dbreg:::projection_memory(c(12500, 1563), 3, length(fixture[["h"]])) < 512 * 1024^2)
  solver = dbreg:::fe_solver(sys, new.env(parent = emptyenv()))
  expect_identical(solver[["method"]], "sparse_cholesky")
  truth = matrix(rnorm(sys[["n"]] * 3L), sys[["n"]], 3L)
  rhs = dbreg:::fe_multiply(sys, truth)
  solution = solver[["solve"]](solver[["project"]](rhs))
  expect_equal(dbreg:::fe_multiply(sys, solution), rhs, tolerance = 1e-10)
  # Independently eliminate the largest FE by pair differencing, then use
  # sparse QR for the remaining incidence design (not the Schur factor).
  big = 12500L
  small = 1563L
  d = data.frame(g = fixture[["i"]], h = fixture[["j"]], w = fixture[["h"]], x = rnorm(2L * big))
  d[["y"]] = 0.3 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(nrow(d), sd = 0.1)
  incidence = Matrix::sparseMatrix(i = rep(seq_len(big), 2L), j = d[["h"]],
    x = c(rep(1, big), rep(-1, big)), dims = c(big, small))
  delta = cbind(x = d[["x"]][seq_len(big)] - d[["x"]][big + seq_len(big)],
    y = d[["y"]][seq_len(big)] - d[["y"]][big + seq_len(big)])
  residual = as.matrix(Matrix::qr.resid(Matrix::qr(incidence[, -1L]), delta))
  beta = sum(residual[, 1L] * residual[, 2L]) / sum(residual[, 1L]^2)
  covariance = sum((residual[, 2L] - beta * residual[, 1L])^2) /
    (nrow(d) - big - small) / sum(residual[, 1L]^2)
  fit = dbreg::dbreg(y ~ x | g + h, data = d, weights = ~w, strategy = "demean")
  expect_equal(unname(stats::coef(fit)), beta, tolerance = 1e-10)
  if (!isTRUE(all.equal(as.numeric(stats::vcov(fit)), covariance, tolerance = 1e-10))) {
    str(list(actual = as.numeric(stats::vcov(fit)), reference = covariance), give.attr = TRUE)
  }
  expect_equal(as.numeric(stats::vcov(fit)), covariance, tolerance = 1e-10)
  expect_identical(fit[["ap_method"]], "bounded_moments")

  # A changing Poisson information matrix reuses the sparse pattern while
  # retaining the same estimates, IID covariance and likelihood as fixest.
  set.seed(993)
  big = 1000L
  small = 600L
  d = data.frame(g = rep(seq_len(big), 4L),
    h = c(rep(seq_len(small), length.out = big), sample.int(small, 3L * big, TRUE)),
    x = rnorm(4L * big), w = runif(4L * big, 0.5, 1))
  d[["y"]] = rpois(nrow(d), exp(0.2 * d[["x"]] + 0.1 * sin(d[["g"]]) + 0.1 * cos(d[["h"]]))) + 1
  fit = dbreg::dbglm(y ~ x | g + h, data = d, weights = ~w, family = poisson(), fixef_rm = "none")
  reference = fixest::feglm(y ~ x | g + h, data = d, weights = ~w, family = poisson(),
    fixef.rm = "none", glm.tol = 1e-12, fixef.tol = 1e-10)
  expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(reference)), tolerance = 1e-9)

})
