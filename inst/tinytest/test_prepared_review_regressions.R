# Prepared structural fallbacks must preserve the same finite-solution and
# offset-start protections as ordinary fits, and release failed refit state.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(id = 1:40, x = rep(0:1, each = 20),
    y = c(rep(c(1, 2, 3, 2), 5), rep(0, 20)))
  DBI::dbWriteTable(conn, "source", d)
  DBI::dbWriteTable(conn, "state", data.frame(id = d[["id"]], w = 1, eta = 0))
  prepared = dbglm_prepare(y ~ x, conn, "source", "id", family = poisson())
  initial = sort(DBI::dbListTables(conn))
  for (warm in c(FALSE, TRUE)) {
    expect_error(dbglm_refit(prepared, "state", etastart = warm), "finite|boundary|separation")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  close(prepared)
  d[["y"]] = c(rep(0, 20), rep(c(0, 1), 10))
  d[["o"]] = rep(c(-100, 0), each = 20)
  DBI::dbWriteTable(conn, "source", d, overwrite = TRUE)
  DBI::dbExecute(conn, "UPDATE state SET eta=CASE WHEN id<=20 THEN -100 ELSE 0 END")
  for (family in list(binomial(), binomial("probit"), binomial("cloglog"))) {
    prepared = dbglm_prepare(y ~ 1, conn, "source", "id", family = family, offset = "o")
    ref = suppressWarnings(glm(y ~ 1, d, offset = o, family = family, start = 0))
    for (warm in c(FALSE, TRUE)) {
      fit = dbglm_refit(prepared, "state", etastart = warm)
      expect_true(fit[["converged"]] && fit[["finite_mle"]])
      expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-7)
      close(fit)
    }
    close(prepared)
    expect_identical(sort(DBI::dbListTables(conn)), c("source", "state"))
  }
})

# dbglm_control() fe_max_iter follows dbreg.cg_max_iter, not dbreg.ap_max_iter.
local({
  saved = options(dbreg.ap_max_iter = 5L, dbreg.cg_max_iter = NULL)
  on.exit(options(saved), add = TRUE)
  expect_identical(dbglm_control()[["fe_max_iter"]], 10000L)
  options(dbreg.cg_max_iter = 250L)
  expect_identical(dbglm_control()[["fe_max_iter"]], 250L)
})
