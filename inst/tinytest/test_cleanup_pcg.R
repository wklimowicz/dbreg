local({
  set.seed(1004)
  n = 2400L
  d = data.frame(g = rep(1:60, each = 40), x = rnorm(n), z = rnorm(n),
    w = sample(c(0.1, 1, 10, 100), n, TRUE))
  d[["h"]] = ((d[["g"]] - 1L) %/% 20L) * 10L + sample.int(10L, n, TRUE)
  d[["y"]] = 0.3 * d[["x"]] - 0.2 * d[["z"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = 0.1)
  d = d[sample.int(n), ]
  reference = fixest::feols(y ~ x + z | g + h, data = d, weights = ~w, vcov = ~h,
    fixef.rm = "none", fixef.tol = 1e-10, ssc = fixest::ssc(K.fixef = "full"))
  conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "pcg_sample", d)
  old = options(dbreg.duckdb_ap_method = "cg", dbreg.glm_fe_cross_max = getOption("dbreg.glm_fe_cross_max"))
  on.exit(options(old), add = TRUE)
  for (cross_limit in c(0, 5e7)) {
    options(dbreg.glm_fe_cross_max = cross_limit)
    fit = dbreg::dbreg(y ~ x + z | g + h, conn = conn, table = "pcg_sample", weights = ~w,
      vcov = ~h, strategy = "demean")
    expect_identical(fit[["ap_method"]], "bounded_pcg")
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
    expect_identical(DBI::dbListTables(conn), "pcg_sample")
  }
  options(dbreg.glm_fe_cross_max = old[["dbreg.glm_fe_cross_max"]])
  shifted = d
  shifted[["x"]] = shifted[["x"]] + 1e9
  shifted[["y"]] = shifted[["y"]] + 2e9
  reference_data = shifted
  reference_data[["x"]] = shifted[["x"]] - 1e9
  reference_data[["y"]] = shifted[["y"]] - 2e9
  shifted_reference = fixest::feols(y ~ x + z | g + h, data = reference_data, weights = ~w, vcov = ~h,
    fixef.rm = "none", fixef.tol = 1e-10, ssc = fixest::ssc(K.fixef = "full"))
  DBI::dbWriteTable(conn, "pcg_shifted", shifted)
  shifted_fit = dbreg::dbreg(y ~ x + z | g + h, conn = conn, table = "pcg_shifted", weights = ~w,
    vcov = ~h, strategy = "demean")
  expect_identical(shifted_fit[["ap_method"]], "bounded_pcg")
  expect_equal(stats::coef(shifted_fit), stats::coef(shifted_reference), tolerance = 1e-8)
  expect_equal(as.numeric(stats::vcov(shifted_fit)), as.numeric(stats::vcov(shifted_reference)), tolerance = 1e-7)
  DBI::dbRemoveTable(conn, "pcg_shifted")
  # A direct projection exposes its preconditioner and component count. This
  # also checks that caller cleanup releases every returned physical relation.
  for (cross_limit in c(0, 5e7)) {
    options(dbreg.glm_fe_cross_max = cross_limit)
    projected = dbreg:::dbreg_duckdb_cg_projections(conn, "FROM pcg_sample", c("g", "h"),
      "y", c("x", "z"), c("x", "z"), "w", scales = c(1, 1, 1), tol = 1e-8)
    expect_equal(projected[["components"]], 3L)
    expect_identical(projected[["preconditioner"]], if (cross_limit == 0) "block_sgs" else "schur_sgs")
    rows = DBI::dbGetQuery(conn, paste("SELECT * FROM", projected[["table"]]))
    expect_equal(nrow(rows), n)
    raw = as.matrix(d[c("y", "x", "z")])
    within = as.matrix(rows[c("y_tilde", "x_tilde", "z_tilde")])
    expect_equal(unname(crossprod(raw, d[["w"]] * raw) - projected[["cross_reduction"]]),
      unname(crossprod(within, rows[["__w"]] * within)), tolerance = 1e-8)
    for (column in c("x_tilde", "z_tilde", "y_tilde")) {
      means = tapply(rows[["__w"]] * rows[[column]], rows[["g"]], sum) / tapply(rows[["__w"]], rows[["g"]], sum)
      expect_true(max(abs(means)) < 1e-8)
    }
    for (table in rev(projected[["all_tables"]])) DBI::dbRemoveTable(conn, table)
    expect_identical(DBI::dbListTables(conn), "pcg_sample")
  }
})
