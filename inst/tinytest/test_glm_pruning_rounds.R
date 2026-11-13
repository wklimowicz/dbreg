local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  make_table = function(sql, materialize) {
    stopifnot(isTRUE(materialize))
    name = basename(tempfile("pruning_test_"))
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", name, "AS", sql))
    name
  }
  drop = function(name) DBI::dbExecute(conn, paste("DROP TABLE", name))
  panel = data.frame(id = 1:7, g = c(1, 1, 2, 2, 1, 3, 3),
    h = c(1, 2, 1, 2, 3, 3, 4), y = c(0, 1, 1, 0, 1, 0, 1))
  DBI::dbWriteTable(conn, "panel", panel)
  owned = make_table("SELECT * FROM panel", TRUE)
  pruned = dbreg:::dbglm_prune(conn, owned, c("g", "h"), "y", "perfect_fit",
    make_table, drop, family = "gaussian")
  expect_equal(pruned[["rounds"]], 3L)
  expect_equal(pruned[["counts"]], c(singleton = 3, infinite_coef = 0))
  expect_equal(DBI::dbReadTable(conn, pruned[["table"]])[["id"]], 1:4)
  expect_equal(sort(DBI::dbListTables(conn)), sort(c("panel", pruned[["table"]])))
  drop(pruned[["table"]])

  panel = data.frame(id = 1:6, g = rep(1:3, each = 2L), y = c(1, 2, 0, 0, 3, 3))
  DBI::dbWriteTable(conn, "panel", panel, overwrite = TRUE)
  owned = make_table("SELECT * FROM panel", TRUE)
  pruned = dbreg:::dbglm_prune(conn, owned, "g", "y", "infinite_coef",
    make_table, drop, family = "poisson")
  expect_equal(pruned[["rounds"]], 1L)
  expect_equal(pruned[["counts"]], c(singleton = 0, infinite_coef = 2))
  expect_equal(sort(DBI::dbReadTable(conn, pruned[["table"]])[["id"]]), c(1L, 2L, 5L, 6L))
})
