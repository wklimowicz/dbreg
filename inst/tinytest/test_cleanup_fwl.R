local({
  set.seed(1006)
  n = 2400L
  d = data.frame(g = rep(1:24, each = 100), h = sample.int(12, n, TRUE), x = rnorm(n),
    w = runif(n, 0.25, 1.5))
  d[["y"]] = 0.4 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = 0.2)
  conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "fwl_source", d)
  reference = fixest::feols(y ~ x | g + h, data = d, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
  # These ordinary IID fits must not request any fit-owned table manager.
  # An accidental return to the old N-row materialization path fails here.
  trace("fit_table_manager", where = asNamespace("dbreg"),
    tracer = quote(stop("Read-only FWL requested a table manager.")), print = FALSE)
  on.exit(untrace("fit_table_manager", where = asNamespace("dbreg")), add = TRUE)
  sources = c("fwl_source",
    "(SELECT * EXCLUDE(g),CAST(9007199254740992 AS BIGINT)+g AS g FROM fwl_source) exact_keys",
    "(SELECT * EXCLUDE(g),1.0+g*2.220446049250313e-16 AS g FROM fwl_source) adjacent_keys",
    "(SELECT * EXCLUDE(g),'key café ' || CAST(g AS VARCHAR) AS g FROM fwl_source) text_keys")
  for (table in sources) {
    fit = dbreg::dbreg(y ~ x | g + h, conn = conn, table = table, strategy = "demean")
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-9)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-8)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
    expect_identical(DBI::dbListTables(conn), "fwl_source")
  }
  untrace("fit_table_manager", where = asNamespace("dbreg"))
  # Retention changes the saved state, while estimation uses the same compact
  # moments. Inference deliberately retains original rows for later requests.
  for (retain in c("prediction", "inference", "all")) {
    fit = dbreg::dbreg(y ~ x | g + h, conn = conn, table = "fwl_source",
      strategy = "demean", retain = retain)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-9)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-8)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
    if (retain %in% c("prediction", "all")) {
      expect_equal(stats::predict(fit, d), as.numeric(stats::predict(reference)), tolerance = 1e-7)
    }
    if (retain %in% c("inference", "all")) {
      expect_equal(as.numeric(stats::vcov(fit, vcov = ~h)),
        as.numeric(stats::vcov(reference, vcov = ~h)), tolerance = 1e-7)
    }
    close(fit)
    expect_identical(DBI::dbListTables(conn), "fwl_source")
  }
  # Uniform rescaling preserves IID and sandwich covariance as well as slopes.
  for (scale in c(1e-200, 1, 1e200)) for (kind in list("iid", "hc1", ~h)) {
    reference = fixest::feols(y ~ x | g + h, data = d, weights = ~w,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), vcov = kind)
    weighted = paste0("(SELECT * EXCLUDE(w),w*(", dbreg:::dbglm_number(scale), ") AS w FROM fwl_source) scaled")
    fit = dbreg::dbreg(y ~ x | g + h, conn = conn, table = weighted, weights = ~w,
      strategy = "demean", vcov = kind)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
    expect_identical(DBI::dbListTables(conn), "fwl_source")
  }
  # Aliased source identifiers must not collide with weights or within columns.
  for (variant in c("w", "__w", "x_tilde")) {
    special = d
    special[["analysis_weight"]] = d[["w"]]
    special[[variant]] = d[["g"]]
    fml = stats::as.formula(paste0("y~x|`", variant, "`+h"))
    DBI::dbWriteTable(conn, "fwl_names", special)
    reference = fixest::feols(fml, data = special, weights = ~analysis_weight,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), vcov = ~h)
    fit = dbreg::dbreg(fml, conn = conn, table = "fwl_names", weights = ~analysis_weight, strategy = "demean", vcov = ~h)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-8)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
    DBI::dbRemoveTable(conn, "fwl_names")
  }
  # Read-only row recomputation handles FWL cancellation without N-row tables.
  d[["near"]] = sin(d[["g"]]) + 1e-5 * d[["x"]]
  d[["absorbed"]] = sin(d[["g"]])
  DBI::dbWriteTable(conn, "fwl_cancellation", d)
  for (fml in list(y ~ near | g + h, y ~ x + absorbed | g + h)) {
    reference = fixest::feols(fml, data = d, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    fit = dbreg::dbreg(fml, conn = conn, table = "fwl_cancellation", strategy = "demean")
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-7)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-6)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
  }
})
