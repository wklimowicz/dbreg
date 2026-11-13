# Separation is a property of outcome support, not a relative tolerance on a
# nearly vanishing score. This one-column finite problem is well conditioned
# after diagonal scaling, even though its positive support is tiny.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  epsilon = 1e-14
  d = data.frame(y = c(rep(1, 20), rep(0, 20)),
    x = c(rep(c(-epsilon, epsilon), 10), rep(1, 20)))
  expected = uniroot(function(b) -epsilon * sinh(epsilon * b) - exp(b),
    c(-100, 0), tol = 1e-13)[["root"]]
  for (strategy in c("moments", "compress")) {
    fit = dbglm(y ~ 0 + x, conn = conn, data = d, family = poisson(), strategy = strategy)
    expect_true(fit[["converged"]] && fit[["finite_mle"]])
    expect_equal(unname(coef(fit)), expected, tolerance = 1e-4)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})

# Every FE group contains positive Poisson counts, yet the zero-only bridge
# can decrease forever while both positive cells remain fixed. A nuisance
# correction must not disappear merely because there are no structural slopes.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(g = rep(c(1, 1, 2), each = 8), h = rep(c(1, 2, 2), each = 8),
    y = rep(c(1, 0, 1), each = 8), x = rep(c(-1, 1), 12), w = 1)
  for (i in 1:2) {
    strategy = c("moments", "compress")[[i]]
    control = dbglm_control(projection_method = c("auto", "sql")[[i]], glm_max_iter = 50L)
    for (fml in list(y ~ 1 | g + h, y ~ x | g + h)) {
      expect_error(dbglm(fml, conn = conn, data = d, family = poisson(), strategy = strategy,
        control = control), "finite|boundary|separation|converge")
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
    # A tiny positive weighted count restores a finite bridge mean. Comparing
    # against its analytic mean also detects premature small-score convergence.
    finite = d
    finite[["y"]][9L] = 1
    finite[["w"]][9L] = 1e-12
    fit = dbglm(y ~ 1 | g + h, conn = conn, data = finite, weights = "w", family = poisson(),
      strategy = strategy, control = control, retain = "inference")
    expected = 1e-12 / (7 + 1e-12)
    expect_true(fit[["converged"]] && fit[["finite_mle"]])
    expect_equal(log(fitted(fit)[9:16]), rep(log(expected), 8), tolerance = 1e-4)
    expect_equal(fitted(fit)[9:16] / expected, rep(1, 8), tolerance = 1e-4)
    close(fit)
    # A weak warm predictor may start after the ordinary projection has already
    # lost this bridge direction. Its independent directional score still gates.
    finite[["eta"]] = ifelse(finite[["g"]] == 1 & finite[["h"]] == 2, -30, 0)
    warm = dbglm(y ~ 1 | g + h, conn = conn, data = finite, weights = "w", family = poisson(),
      etastart = "eta", strategy = strategy, control = control, retain = "inference")
    expect_true(warm[["converged"]] && warm[["finite_mle"]])
    expect_equal(fitted(warm)[9:16] / expected, rep(1, 8), tolerance = 1e-4)
    close(warm)
    d[["eta"]] = finite[["eta"]]
    expect_error(dbglm(y ~ 1 | g + h, conn = conn, data = d, family = poisson(),
      etastart = "eta", strategy = strategy, control = control), "finite|boundary|separation|converge")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})

# Two mixed-outcome blocks joined by one zero-only bridge have a binomial
# recession direction although recursive pure-group pruning removes no rows.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  cells = data.frame(g = c(1, 1, 2, 2, 3, 3, 4, 4, 1),
    h = c(1, 2, 1, 2, 3, 4, 3, 4, 3), y = c(0, 1, 1, 0, 0, 1, 1, 0, 0))
  d = cells[rep(seq_len(nrow(cells)), each = 5), ]
  families = list(binomial(), binomial("probit"), binomial("cloglog"))
  for (i in seq_along(families)) {
    strategy = if (i == 1L) "moments" else "compress"
    control = dbglm_control(projection_method = if (i == 2L) "sql" else "auto", glm_max_iter = 50L)
    expect_error(dbglm(y ~ 1 | g + h, conn = conn, data = d, family = families[[i]],
      strategy = strategy, control = control), "finite|boundary|separation|converge")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
    if (i == 3L) {
      positive = d
      positive[["y"]][(nrow(positive) - 4L):nrow(positive)] = 1
      positive[["eta"]] = c(rep(0, nrow(positive) - 5L), rep(4, 5))
      expect_error(dbglm(y ~ 1 | g + h, conn = conn, data = positive,
        family = families[[i]], etastart = "eta", strategy = strategy, control = control),
        "finite|boundary|separation|converge")
      expect_identical(sort(DBI::dbListTables(conn)), initial)
    }
    finite = d
    finite[["y"]][nrow(finite)] = 1
    fit = dbglm(y ~ 1 | g + h, conn = conn, data = finite, family = families[[i]],
      strategy = strategy, control = control, retain = "inference")
    expect_true(fit[["converged"]] && fit[["finite_mle"]])
    expect_equal(tail(fitted(fit), 5), rep(.2, 5), tolerance = 1e-8)
    close(fit)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
})
