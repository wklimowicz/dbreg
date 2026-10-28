# Initial overflow is an invalid information matrix, never collinearity.
local({
  for (H in list(matrix(NaN, 1L, 1L), matrix(c(1, Inf, Inf, 1), 2L))) {
    expect_error(dbreg:::dbglm_basis(H, 1e-10), "initial information is nonfinite")
  }
  expect_identical(dbreg:::dbglm_basis(matrix(c(1, 0, 0, 0), 2L), 1e-10), 1L)
  expect_identical(dbreg:::dbglm_basis(matrix(1, 2L, 2L), 1e-10), 1L)
  expect_error(dbreg:::dbglm_basis(matrix(0, 0L, 0L), 1e-10), "No identified")
  expect_error(dbreg:::dbglm_basis(matrix(0, 2L, 2L), 1e-10), "No identified")
})

local({
  set.seed(28415)
  n = 120L
  z = stats::rnorm(n, sd = 0.4)
  y = 0.3 + 0.7 * z + stats::rnorm(n, sd = 0.2)
  X = cbind(1, z)
  inverse = solve(crossprod(X))
  beta = inverse %*% crossprod(X, y)
  residual = y - drop(X %*% beta)
  V = inverse * sum(residual^2) / (n - ncol(X))
  target = drop(X %*% beta)
  # Disable the optional K/G adjustments for a plain lm covariance oracle.
  plain = list(K.adj = FALSE, G.adj = FALSE)
  for (scale in c(1, 1e100)) {
    d = data.frame(x = z * scale, y = y)
    fit = dbreg::dbglm(y ~ x, data = d, ssc = plain)
    units = c(1, scale)
    expect_equal(as.numeric(stats::coef(fit)) * units, as.numeric(beta), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit) * outer(units, units)), as.numeric(V), tolerance = 1e-8)
    expect_equal(stats::predict(fit, d), target, tolerance = 1e-8)
    expect_identical(fit[["collin.var"]], character())
  }
  finite = data.frame(x = z, copy = z, zero = 0, y = y)
  fit = dbreg::dbglm(y ~ x + copy + zero, data = finite, ssc = plain)
  expect_equal(as.numeric(stats::coef(fit)), as.numeric(beta), tolerance = 1e-8)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(V), tolerance = 1e-8)
  expect_true(all(is.na(stats::coef(fit, collin = TRUE)[c("copy", "zero")])))
  expect_equal(stats::predict(fit, finite), target, tolerance = 1e-8)

  huge = data.frame(x = z * 1e200, y = y)
  # Capture the real owned connection at cleanup entry, then inspect its lifetime.
  capture = new.env(parent = emptyenv())
  assign("dbreg_initial_guard_capture__", capture, envir = .GlobalEnv)
  on.exit(rm("dbreg_initial_guard_capture__", envir = .GlobalEnv), add = TRUE)
  trace("release_fit_resources", where = asNamespace("dbreg"), print = FALSE,
    tracer = quote({
      state = get("dbreg_initial_guard_capture__", envir = .GlobalEnv)
      if (isTRUE(resources[["own_conn"]])) state[["conn"]] = resources[["conn"]]
    }))
  on.exit(untrace("release_fit_resources", where = asNamespace("dbreg")), add = TRUE)
  expect_error(dbreg::dbglm(y ~ x, data = huge), "initial information is nonfinite")
  expect_true(!is.null(capture[["conn"]]))
  expect_true(!DBI::dbIsValid(capture[["conn"]]))

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "caller_sample", huge)
  expect_error(dbreg::dbglm(y ~ x, conn = conn, table = "caller_sample"), "initial information is nonfinite")
  expect_true(DBI::dbIsValid(conn))
  expect_identical(DBI::dbListTables(conn), "caller_sample")
  expect_equal(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM caller_sample")[["n"]], n)
})
