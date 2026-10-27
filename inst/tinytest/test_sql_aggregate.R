# Ordinary generated sums retain their own policy, separate from working sums.
library(DBI)
source(system.file("tinytest/helper_glm_databricks.R", package = "dbreg"), local = TRUE)
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", temp_directory = "")))
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit({
    dbreg:::release_fit_resources(resources)
    DBI::dbDisconnect(conn, shutdown = TRUE)
  })
  remote = dbglm_test_databricks_connection(conn)
  expect_identical(dbreg:::sql_sum(conn, c("w", "w*(x)")), c("FSUM(w)", "FSUM(w*(x))"))
  expect_identical(dbreg:::sql_sum(remote, c("w", "w*(x)")), c("SUM(w)", "SUM(w*(x))"))
  expect_identical(dbreg:::sql_sum(conn, character()), character())
  for (backend in c("mysql", "sqlserver", "bigquery", "athena", "other")) {
    expect_identical(dbreg:::sql_sum(NULL, c("x", "x*y"), list(name = backend)), c("SUM(x)", "SUM(x*y)"))
    expect_identical(dbreg:::sql_stable_sum(NULL, "x", list(name = backend)), "SUM(x)")
  }
  expect_false(grepl("DECIMAL", dbreg:::sql_sum(remote, "x"), fixed = TRUE))
  expect_true(grepl("DECIMAL(38,18)", dbreg:::sql_stable_sum(remote, "x"), fixed = TRUE))
  expect_false(exists("backend_sql", envir = asNamespace("dbreg"), inherits = FALSE))
  # Exact production SQL shapes, including aliases/order, preserve DuckDB plans.
  expect_identical(dbreg:::fe_value_terms(conn, c("x", "y")),
    c("FSUM(w) AS __h", "FSUM(w*(x)) AS __b1", "FSUM(w*(y)) AS __b2"))
  expect_identical(dbreg:::fe_value_terms(remote, c("x", "y")),
    c("SUM(w) AS __h", "SUM(w*(x)) AS __b1", "SUM(w*(y)) AS __b2"))
  expect_identical(dbreg:::fwl_information_terms(conn, c("y", "x")), c(
    "COUNT(*) AS fwl_n", "FSUM(w) AS fwl_sw", "FSUM(w*(y)*(y)) AS h_1_1",
    "FSUM(w*(x)*(y)) AS h_2_1", "FSUM(w*(x)*(x)) AS h_2_2"))
  # Independent binary-exact sum: two unit terms survive large cancellation.
  cancellation = DBI::dbGetQuery(conn, paste0("SELECT ", dbreg:::sql_sum(conn, "v"),
    " AS value FROM (VALUES (1e16),(1.0),(1.0),(-1e16)) t(v)"))
  expect_equal(cancellation[["value"]], 2, tolerance = 0)
  for (engine in list(conn, remote)) {
    mixed = DBI::dbGetQuery(engine, paste0("SELECT ", dbreg:::sql_sum(engine, "v"),
      " AS value FROM (VALUES (1.0),(NULL),(2.0)) t(v)"))
    empty = DBI::dbGetQuery(engine, paste0("SELECT ", dbreg:::sql_sum(engine, "v"),
      " AS value FROM (VALUES (1.0)) t(v) WHERE FALSE"))
    expect_equal(mixed[["value"]], 3)
    expect_true(is.na(empty[["value"]]))
  }

  set.seed(8631)
  d = expand.grid(g = seq_len(3L), h = seq_len(2L), replicate = seq_len(4L))
  d[["w"]] = sample(1:3, nrow(d), TRUE)
  d[["x"]] = stats::rnorm(nrow(d))
  d[["y"]] = .3*d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + stats::rnorm(nrow(d))
  DBI::dbWriteTable(conn, "aggregate_fixture", d)
  # Actual compact FE summary producer vs independent weighted dummy products.
  dummy = cbind(stats::model.matrix(~factor(g)-1, d), stats::model.matrix(~factor(h)-1, d))
  for (engine in list(conn, remote)) {
    moments = dbreg:::query_fe_moments(engine, "aggregate_fixture", c("g", "h"), c("x", "y"))
    expect_equal(unname(moments[["B"]]), unname(crossprod(dummy, d[["w"]] * as.matrix(d[c("x", "y")]))), tolerance = 1e-12)
    expect_equal(moments[["sys"]][["diag"]], as.numeric(crossprod(dummy, d[["w"]])), tolerance = 1e-12)
  }
  # General weighted FWL matrix independently recovered with weighted QR.
  nuisance = stats::model.matrix(~factor(g)+factor(h), d)
  within = qr.resid(qr(nuisance * sqrt(d[["w"]])), as.matrix(d[c("y", "x")]) * sqrt(d[["w"]]))
  reference_M = crossprod(within)
  reference_beta = reference_M["x", "y"] / reference_M["x", "x"]
  reference_rss = sum((within[, "y"] - reference_beta * within[, "x"])^2)
  fit = dbreg::dbreg(y ~ x | g + h, conn = conn, table = "aggregate_fixture", weights = "w", strategy = "demean")
  expect_equal(unname(stats::coef(fit)), reference_beta, tolerance = 1e-10)
  expect_equal(attr(stats::vcov(fit), "rss"), reference_rss, tolerance = 1e-10)
  reference = stats::lm(y ~ x + factor(g) + factor(h), d, weights = w)
  expect_equal(as.numeric(stats::vcov(fit)), as.numeric(stats::vcov(reference)["x", "x"]), tolerance = 1e-10)
  close(fit)
  # Actual two-FE recovery means/certificate sums, not a detached toy SELECT.
  manager = dbreg:::fit_table_manager(conn, resources)
  projection = dbreg:::prediction_two_fe_projection(conn, "aggregate_fixture", c("g", "h"), c("y", "x"),
    manager[["create"]], manager[["drop"]])
  projected = DBI::dbGetQuery(conn, paste("SELECT g,h,w,r1,r2 FROM", projection[["residual"]]))
  membership = cbind(stats::model.matrix(~factor(g)-1, projected), stats::model.matrix(~factor(h)-1, projected))
  expect_true(max(abs(crossprod(membership, projected[["w"]] * as.matrix(projected[c("r1", "r2")])))) < 1e-8)
  expect_equal(crossprod(as.matrix(projected[c("r1", "r2")]) * sqrt(projected[["w"]])) *
    outer(projection[["scales"]], projection[["scales"]]), reference_M, tolerance = 1e-8, check.attributes = FALSE)
})
