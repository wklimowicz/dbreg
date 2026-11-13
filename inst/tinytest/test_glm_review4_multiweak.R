# Two independent weak bridges defeat a single saved-direction certificate:
# deliberately choose opposite warm errors whose directional scores cancel.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  cells = data.frame(g = c(1, 2, 3, 1, 2), h = c(1, 2, 3, 2, 3), y = c(1, 1, 1, 0, 0))
  d = cells[rep(seq_len(nrow(cells)), each = 20), ]
  d[["w"]] = 1
  d[["y"]][c(61L, 81L)] = 1
  d[["w"]][c(61L, 81L)] = 1e-12
  expected = 1e-12 / (19 + 1e-12)
  alpha = log(expected)
  error1 = log(1.5)
  error2 = uniroot(function(value) (alpha + error1) * (1 - exp(error1)) +
    (alpha + value) * (1 - exp(value)), c(-5, -1e-4), tol = 1e-13)[["root"]]
  d[["eta"]] = c(rep(0, 60), rep(alpha + error1, 20), rep(alpha + error2, 20))
  for (strategy in c("moments", "compress")) for (method in c("auto", "sql")) {
    for (warm in c(FALSE, TRUE)) {
      fit = dbglm(y ~ 1 | g + h, conn = conn, data = d, weights = "w", family = poisson(),
        etastart = if (warm) "eta", strategy = strategy, retain = "inference",
        control = dbglm_control(projection_method = method))
      expect_true(fit[["converged"]] && fit[["finite_mle"]])
      expect_equal(fitted(fit)[c(61L, 81L)] / expected, c(1, 1), tolerance = 1e-6)
      expect_equal(fitted(fit)[1:60], rep(1, 60), tolerance = 1e-10)
      expect_identical(fit[["projection_method"]], "sql_qr")
      close(fit)
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
  }
})

# Numerical recovery must obey the existing transfer plan and clean up on a
# capacity error; a weak fit may not silently fall back to a zero correction.
# Without room for QR, the CG fallback meets the same transfer budget.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(g = rep(c(1, 1, 2), each = 8), h = rep(c(1, 2, 2), each = 8),
    y = rep(c(1, 0, 1), each = 8), w = 1)
  d[["y"]][9L] = 1
  d[["w"]][9L] = 1e-12
  previous = options(dbreg.projection_max_rows = 3L)
  on.exit(options(previous), add = TRUE)
  expect_error(dbglm(y ~ 1 | g + h, conn = conn, data = d, weights = "w", family = poisson(),
    strategy = "moments", control = dbglm_control(projection_method = "sql")), "transfer budget")
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
