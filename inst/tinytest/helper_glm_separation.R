# Known separating directions give an oracle independent of either fitter.
# A successful finite-MLE result is never acceptable on these designs.
dbglm_separation_cases = function() {
  complete = data.frame(x = rep(c(-2, -1, 1, 2), each = 4), y = rep(c(0, 0, 1, 1), each = 4))
  quasi = expand.grid(x = c(-1, 0, 1), replicate = 0:1, f = 1:3, g = 1:2)
  quasi[["y"]] = ifelse(quasi[["x"]] == 0, quasi[["replicate"]], as.integer(quasi[["x"]] > 0))
  combined = expand.grid(x = c(-4, -1, 1, 4), z = c(-2, 2), f = 1:3, g = 1:2)
  direction = 2 * combined[["x"]] + combined[["z"]] + combined[["f"]] - combined[["g"]] + .25
  combined[["y"]] = as.integer(direction > 0)
  list(
    complete = list(data = complete, formula = y ~ x, direction = complete[["x"]]),
    quasi = list(data = quasi, formula = y ~ x, direction = quasi[["x"]]),
    quasi_fe = list(data = quasi, formula = y ~ x | f + g, direction = quasi[["x"]]),
    combined_fe = list(data = combined, formula = y ~ x + z | f + g, direction = direction))
}

dbglm_test_separation = function(conn, verbose = FALSE) {
  cases = dbglm_separation_cases()
  for (name in names(cases)) {
    if (verbose) message("[GLM separation] ", name)
    case = cases[[name]]
    d = case[["data"]]
    margin = (2 * d[["y"]] - 1) * case[["direction"]]
    expect_true(all(margin >= 0) && any(margin > 0), info = name)
    if ("f" %in% names(d)) for (f in c("f", "g")) {
      expect_true(all(vapply(split(d[["y"]], d[[f]]), function(y) length(unique(y)) == 2L, logical(1))), info = name)
    }
    before = dbglm_validation_tables(conn)
    table = dbglm_validation_upload(conn, d)
    export = basename(tempfile("glm_separated_sample_"))
    tryCatch({
      expect_error(dbglm(case[["formula"]], conn = conn, table = table, family = binomial("cloglog"),
        fixef_rm = "none", sample_table = export),
        "boundary|converge|singular|halving|separation|identification", info = name)
      expect_identical(dbglm_validation_tables(conn), sort(c(before, table)), info = name)
    }, finally = {dbglm_validation_drop(conn, export); dbglm_validation_drop(conn, table)})
    expect_identical(dbglm_validation_tables(conn), before, info = name)
  }

  # One contrary outcome at each covariate value makes the MLE finite. The
  # saturated two-point solution is analytic, including its information matrix.
  finite = expand.grid(x = c(-1, 1), y = 0:1)
  finite[["a"]] = ifelse(finite[["y"]] == as.integer(finite[["x"]] > 0), 1000L, 1L)
  probabilities = c(1, 1000) / 1001
  family = binomial("cloglog")
  eta = family[["linkfun"]](probabilities)
  expected = c("(Intercept)" = mean(eta), x = diff(eta) / 2)
  D = cbind(1, c(-1, 1))
  H = crossprod(D, (1001 * family[["mu.eta"]](eta)^2 / (probabilities * (1 - probabilities))) * D)
  before = dbglm_validation_tables(conn)
  table = dbglm_validation_upload(conn, finite)
  tryCatch({
    fit = dbglm(y ~ x, conn = conn, table = table, family = family, weights = "a", ssc = list(K.adj = FALSE))
    expect_equal(coef(fit), expected, tolerance = 1e-7)
    expect_equal(unname(vcov(fit)), unname(solve(H)), tolerance = 1e-7)
    expect_true(fit[["converged"]] && fit[["score_max"]] <= fit[["control"]][["score_tol"]])
    # The same finite slope survives absorption, with both outcomes in every FE.
    crossed = merge(finite, expand.grid(f = 1:3, g = 1:2))
    absorbed = dbglm_validation_upload(conn, crossed)
    tryCatch({
      fit_fe = dbglm(y ~ x | f + g, conn = conn, table = absorbed, weights = "a", ssc = list(K.adj = FALSE))
      expect_equal(coef(fit_fe)[["x"]], expected[["x"]], tolerance = 1e-7)
      expect_equal(vcov(fit_fe)[1, 1], solve(H)[2, 2] / 6, tolerance = 1e-7)
      expect_true(fit_fe[["converged"]] && fit_fe[["score_max"]] <= fit_fe[["control"]][["score_tol"]])
    }, finally = dbglm_validation_drop(conn, absorbed))
  }, finally = dbglm_validation_drop(conn, table))
  expect_identical(dbglm_validation_tables(conn), before)
}
