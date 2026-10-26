library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  set.seed(1014)
  n = 1800L
  d = data.frame(g = rep(1:18, each = 100), h = sample(1:9, n, TRUE),
    x = stats::rnorm(n), w = stats::runif(n, 0.2, 1.5))
  d[["y"]] = 0.4 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + stats::rnorm(n, sd = 0.2)
  native = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  DBI::dbWriteTable(native, "fwl_native_source", d)
  conn = dbglm_test_databricks_connection(native)
  for (kind in list("iid", "hc1", ~h)) for (retain in c("none", "prediction")) {
    reference = fixest::feols(y ~ x | g + h, data = d, weights = ~w,
      vcov = kind, fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    fit = dbreg(y ~ x | g + h, conn = conn, table = "fwl_native_source", weights = ~w,
      vcov = kind, strategy = "demean", retain = retain)
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-7)
    expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)), tolerance = 1e-7)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
    if (retain == "prediction") expect_equal(stats::predict(fit, d),
      as.numeric(stats::predict(reference)), tolerance = 1e-7)
    expect_identical(DBI::dbListTables(native), "fwl_native_source")
  }
  for (variant in c("text", "logical", "date", "float")) {
    special = d
    special[["g"]] = switch(variant,
      text = paste0("NaN café '", d[["g"]]), logical = d[["g"]] %% 2L == 0L,
      date = as.Date("2026-01-01") + d[["g"]], float = d[["g"]] / 10)
    special[["y"]] = d[["y"]]
    DBI::dbWriteTable(native, "fwl_typed", special, overwrite = TRUE)
    table = if (variant == "float") "(SELECT * EXCLUDE(g),CAST(g AS FLOAT) AS g FROM fwl_typed) float_keys" else "fwl_typed"
    reference = fixest::feols(y ~ x | g + h, data = special, weights = ~w,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    fit = dbreg(y ~ x | g + h, conn = conn, table = table, weights = ~w,
      strategy = "demean", retain = "prediction")
    expect_equal(stats::coef(fit), stats::coef(reference), tolerance = 1e-7)
    expect_equal(stats::predict(fit, special), as.numeric(stats::predict(reference)), tolerance = 1e-7)
    expect_equal(fit[["projection_diagnostics"]][[1L]][["details"]][["observation_tables"]], 0L)
  }
})
