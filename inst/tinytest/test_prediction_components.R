local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  resources = dbreg:::dbglm_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::dbglm_release_resources(resources), add = TRUE)
  manager = dbreg:::dbglm_table_manager(conn, resources)
  DBI::dbWriteTable(conn, "caller_graph", data.frame(g = c(1, 3, 3, 7, 7, 9), h = c(2, 2, 5, 11, 21, 21)))
  labels = dbreg:::prediction_two_fe_components(conn, "caller_graph", c("g", "h"), manager[["create"]], manager[["drop"]])
  left = DBI::dbGetQuery(conn, paste("SELECT * FROM", labels[["left"]], "ORDER BY id"))
  expect_equal(left[["id"]], c(1, 3, 7, 9))
  expect_equal(left[["component"]], c(1, 1, 7, 7))
  # The graph has large diameter, but pointer compression avoids propagating
  # just one edge per iteration. Check the final component membership itself.
  DBI::dbExecute(conn, paste0("CREATE TABLE caller_chain AS SELECT i AS g,i AS h FROM range(1,1001) t(i) ",
    "UNION ALL SELECT i+1 AS g,i AS h FROM range(1,1001) t(i)"))
  labels = dbreg:::prediction_two_fe_components(conn, "caller_chain", c("g", "h"), manager[["create"]], manager[["drop"]])
  values = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,MIN(component) AS lo,MAX(component) AS hi FROM", labels[["left"]]))
  expect_equal(as.numeric(values[["n"]]), 1001)
  expect_equal(as.numeric(values[["lo"]]), 1)
  expect_equal(as.numeric(values[["hi"]]), 1)
  expect_error(dbreg:::prediction_two_fe_components(conn, "caller_chain", c("g", "h"), manager[["create"]],
    manager[["drop"]], max_iter = 1L), "labels did not converge")
  dbreg:::dbglm_release_resources(resources)
  expect_equal(sort(DBI::dbListTables(conn)), c("caller_chain", "caller_graph"))
})
