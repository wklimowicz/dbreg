library(DBI)
library(dbreg)
library(fixest)

## ---- backend capabilities -------------------------------------------------

odbc_info = dbreg:::backend_info_from_metadata(
  info = list(
    dbms.name = "Spark SQL",
    drivername = "Databricks ODBC Driver"
  ),
  connection_string = list(useNativeQuery = "1")
)
expect_equal(odbc_info[["name"]], "spark")
expect_true(odbc_info[["is_databricks"]])
expect_true(odbc_info[["supports_physical_temp_tables"]])
expect_identical(odbc_info[["native_query"]], TRUE)

jdbc_info = dbreg:::backend_info_from_metadata(
  info = list(dbms.name = "Databricks"),
  connection_string = "jdbc:databricks://example;UseNativeQuery=0"
)
expect_true(jdbc_info[["is_databricks"]])
expect_identical(jdbc_info[["native_query"]], FALSE)
expect_error(
  dbreg:::validate_databricks_ap_backend(jdbc_info),
  "useNativeQuery = TRUE"
)

adbc_info = dbreg:::backend_info_from_metadata(
  info = list(dbms.name = "Databricks SQL")
)
expect_true(adbc_info[["is_databricks"]])
expect_true(is.na(adbc_info[["native_query"]]))
expect_silent(dbreg:::validate_databricks_ap_backend(adbc_info))

spark_info = dbreg:::backend_info_from_metadata(
  info = list(dbms.name = "Apache Spark SQL")
)
expect_false(spark_info[["is_databricks"]])
expect_false(spark_info[["supports_physical_temp_tables"]])
expect_error(
  dbreg:::validate_databricks_ap_backend(spark_info),
  "requires Databricks"
)

## ---- bounded window AP ----------------------------------------------------

set.seed(2468)
n = 600L
ap_dat = data.frame(
  fe1 = sample(1:45, n, replace = TRUE),
  fe2 = sample(1:31, n, replace = TRUE),
  fe3 = sample(1:13, n, replace = TRUE),
  cluster = sample(1:20, n, replace = TRUE),
  x1 = rnorm(n),
  x2 = rnorm(n),
  obs_weight = runif(n, 0.2, 2.5)
)
ap_dat = ap_dat[!(ap_dat[["fe1"]] %% 7 == 0 & ap_dat[["fe2"]] %% 5 == 0), ]
ap_dat[["y"]] = 0.7 * ap_dat[["x1"]] - 0.25 * ap_dat[["x2"]] +
  rnorm(45)[ap_dat[["fe1"]]] + rnorm(31)[ap_dat[["fe2"]]] +
  rnorm(13)[ap_dat[["fe3"]]] + rnorm(nrow(ap_dat), sd = 0.4)

con = dbConnect(duckdb::duckdb())
dbWriteTable(con, "databricks_ap_fixture", ap_dat)

run_window_ap = function(chunk_cycles, materialize_source = FALSE) {
  dbreg:::dbreg_window_alternating_projections(
    conn = con,
    from_statement = "FROM databricks_ap_fixture",
    fe = c("fe1", "fe2", "fe3"),
    yvar = "y",
    xvars_sql = c("x1", "x2"),
    xvar_names = c("x1", "x2"),
    weights = "obs_weight",
    cluster_var = "cluster",
    max_iter = 100L,
    tol = 1e-10,
    chunk_cycles = chunk_cycles,
    materialize_source = materialize_source
  )
}

window_fit_1 = run_window_ap(1L)
window_fit_3 = run_window_ap(3L)
window_fit_materialized = run_window_ap(3L, materialize_source = TRUE)

window_coef = function(ap_result) {
  moments = dbGetQuery(con, paste0(
    "SELECT ",
    "SUM(__w * x1_tilde * x1_tilde) AS x11, ",
    "SUM(__w * x1_tilde * x2_tilde) AS x12, ",
    "SUM(__w * x2_tilde * x2_tilde) AS x22, ",
    "SUM(__w * x1_tilde * y_tilde) AS x1y, ",
    "SUM(__w * x2_tilde * y_tilde) AS x2y ",
    "FROM ", ap_result[["table"]]
  ))
  solve(
    matrix(c(
      moments[["x11"]], moments[["x12"]],
      moments[["x12"]], moments[["x22"]]
    ), 2, 2),
    c(moments[["x1y"]], moments[["x2y"]])
  )
}

coef_1 = window_coef(window_fit_1)
coef_3 = window_coef(window_fit_3)
coef_materialized = window_coef(window_fit_materialized)
reference = feols(
  y ~ x1 + x2 | fe1 + fe2 + fe3,
  data = ap_dat,
  weights = ~obs_weight,
  fixef.rm = "none"
)

expect_equal(coef_1, as.numeric(coef(reference)), tolerance = 1e-8)
expect_equal(coef_3, as.numeric(coef(reference)), tolerance = 1e-8)
expect_equal(coef_1, coef_3, tolerance = 1e-10)
expect_equal(coef_1, coef_materialized, tolerance = 1e-10)
expect_true(window_fit_1[["iterations"]] > 3L)
expect_true(window_fit_1[["max_abs"]] < 1e-10)
expect_true(window_fit_3[["iterations"]] %% 3L == 0L)
expect_true(window_fit_3[["max_abs"]] < 1e-10)
expect_equal(length(window_fit_3[["all_tables"]]), 2L)
expect_equal(length(window_fit_materialized[["all_tables"]]), 3L)

for (result in list(window_fit_1, window_fit_3, window_fit_materialized)) {
  for (table_name in result[["all_tables"]]) {
    dbreg:::drop_ap_table_if_exists(con, table_name, result[["backend_info"]])
  }
}
expect_false(any(grepl("^dbreg_ap_", dbListTables(con))))

## ---- validation and failure cleanup --------------------------------------

expect_error(
  dbreg:::dbreg_window_alternating_projections(
    con, "FROM databricks_ap_fixture", c("fe1", "fe2"),
    "y", "x1", "x1", "obs_weight", chunk_cycles = 0L
  ),
  "positive integer"
)
expect_error(
  dbreg:::dbreg_window_alternating_projections(
    con, "FROM databricks_ap_fixture", c("fe1", "fe2"),
    "y", "x1", "x1", "obs_weight", tol = Inf
  ),
  "finite, positive"
)

expect_error(
  dbreg:::dbreg_window_alternating_projections(
    conn = con,
    from_statement = "FROM databricks_ap_fixture",
    fe = c("fe1", "fe2", "fe3"),
    yvar = "y",
    xvars_sql = c("x1", "x2"),
    xvar_names = c("x1", "x2"),
    weights = "obs_weight",
    max_iter = 2L,
    tol = 1e-30,
    chunk_cycles = 3L
  ),
  "did not converge within 2 iterations"
)
expect_false(any(grepl("^dbreg_ap_", dbListTables(con))))

## ---- public diagnostics and cleanup --------------------------------------

public_fit = dbreg(
  y ~ x1 + x2 | fe1 + fe2 + fe3,
  conn = con,
  table = "databricks_ap_fixture",
  weights = "obs_weight",
  strategy = "demean",
  vcov = "hc1"
)
expect_true(public_fit[["ap_iterations"]] >= 1L)
expect_true(public_fit[["ap_max_abs"]] < getOption("dbreg.ap_tol", 1e-10))
expect_false(any(grepl("^dbreg_ap_", dbListTables(con))))

analytic_fit = dbreg(
  y ~ x1 + x2 | fe1,
  conn = con,
  table = "databricks_ap_fixture",
  strategy = "demean"
)
expect_true(is.null(analytic_fit[["ap_iterations"]]))
expect_true(is.null(analytic_fit[["ap_max_abs"]]))

dbDisconnect(con, shutdown = TRUE)
