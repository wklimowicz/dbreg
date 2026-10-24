# Real database resources must survive or be released according to ownership,
# including failures between registering coefficients and materializing them.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "caller_data", data.frame(x = 1:3))

  resources = dbreg:::dbglm_resources()
  resources[["conn"]] = conn
  manager = dbreg:::dbglm_table_manager(conn, resources)
  create = manager[["create"]]

  # Integer64 IDs must be decoded before combining groups, including when
  # several projected columns share one coefficient lookup.
  info = list(list(id = bit64::as.integer64(c(1, 3))), list(id = bit64::as.integer64(2)))
  coefficients = matrix(c(0.2, -0.4, 0.6, 1.2, 1.4, -1.6), ncol = 2)
  frame = dbreg:::dbglm_coefficient_frame(info, coefficients)
  expected = data.frame(op = c(1L, 1L, 2L), group = c(1, 3, 2),
    c1 = c(0.2, -0.4, 0.6), c2 = c(1.2, 1.4, -1.6))
  names(expected)[1:2] = c("__glm_op", "__glm_group")
  expect_equal(frame, expected)
  dictionary = dbreg:::dbglm_coefficient_lookup(conn, frame, create)
  stored = DBI::dbGetQuery(conn, paste("SELECT * FROM", dictionary[["source"]],
    "ORDER BY __glm_op, __glm_group"))
  expect_equal(stored, expected)
  tables_before_failure = sort(DBI::dbListTables(conn))

  fail_create = function(...) stop("injected materialization failure")
  expect_error(dbreg:::dbglm_coefficient_lookup(conn, frame, fail_create),
    "injected materialization failure")
  expect_equal(sort(DBI::dbListTables(conn)), tables_before_failure)

  # A failed CREATE is tracked before execution; cleanup tolerates its absence.
  expect_error(create("SELECT missing_column FROM caller_data"), "missing_column")
  dbreg:::dbglm_release_resources(resources)
  expect_true(DBI::dbIsValid(conn))
  expect_equal(DBI::dbListTables(conn), "caller_data")
  expect_equal(DBI::dbReadTable(conn, "caller_data"), data.frame(x = 1:3))
  dbreg:::dbglm_release_resources(resources)
  expect_equal(DBI::dbListTables(conn), "caller_data")
})
