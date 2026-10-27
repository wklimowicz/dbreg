# Type/exclusion emitters compose actual SQL without inspecting whole statements.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  native = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", temp_directory = "")))
  on.exit(DBI::dbDisconnect(native, shutdown = TRUE))
  remote = dbglm_test_databricks_connection(native)
  expect_identical(dbreg:::sql_text_type(native), "VARCHAR")
  expect_identical(dbreg:::sql_text_type(remote), "STRING")
  expect_identical(dbreg:::sql_cast_text(native, c("NULL", "x")), c("CAST(NULL AS VARCHAR)", "CAST(x AS VARCHAR)"))
  expect_identical(dbreg:::sql_cast_text(remote, c("NULL", "x")), c("CAST(NULL AS STRING)", "CAST(x AS STRING)"))
  expect_identical(dbreg:::sql_cast_text(native, character()), character())
  expect_identical(dbreg:::sql_select_except(native, character()), "*")
  expect_error(dbreg:::sql_select_except(NULL, "x", list(name = "mysql")), "DuckDB and Spark only")
  expect_true(startsWith(dbreg:::sql_select_except(native, c("x", "odd name")), "* EXCLUDE("))
  expect_true(startsWith(dbreg:::sql_select_except(remote, c("x", "odd name")), "* EXCEPT("))
  for (conn in list(native, remote)) {
    query = paste0("SELECT ", dbreg:::sql_select_except(conn, "odd name"), ",",
      dbreg:::sql_cast_text(conn, "CAST(9007199254740993 AS BIGINT)"), " AS exact_key,",
      dbreg:::sql_cast_text(conn, "NULL"), ' AS missing FROM (SELECT 1 AS x,2 AS "odd name",3 AS y) t')
    value = dbreg:::backend_get_query(conn, query)
    expect_identical(names(value), c("x", "y", "exact_key", "missing"))
    expect_identical(value[["exact_key"]], "9007199254740993")
    expect_true(is.na(value[["missing"]]))
    # Native family CTE builders must own preparation-column removal too.
    source = "(SELECT 0.0 AS __dbglm_o,0.4 AS __dbglm_y,1.0 AS __dbglm_a,0.2 AS x) t"
    for (link in c("probit", "cloglog")) {
      cte = dbreg:::dbglm_working_sql(source, "x", 1, link, conn = conn)
      rows = dbreg:::backend_get_query(conn, paste0(cte, "SELECT * FROM working"))
      expect_false(any(c("__dbglm_abs_eta", "__dbglm_eh") %in% names(rows)))
      family = stats::binomial(link)
      expect_equal(rows[["mu"]], family[["linkinv"]](.2), tolerance = 1e-12)
      expect_equal(rows[["d"]], family[["mu.eta"]](.2), tolerance = 1e-12)
    }
  }
})
