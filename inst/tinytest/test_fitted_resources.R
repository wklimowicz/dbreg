local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(if (DBI::dbIsValid(conn)) DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "caller_data", data.frame(x = 1:3))
  DBI::dbWriteTable(conn, "caller_export", data.frame(id = 1:3))
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  manager = dbreg:::fit_table_manager(conn, resources)
  first = manager[["create"]]("SELECT x*2 AS learned FROM caller_data")
  second = manager[["create"]](paste("SELECT * FROM", first))
  unused = manager[["create"]]("SELECT 0 AS intermediate")
  names_owned = resources[["tables"]]
  expect_error(dbreg:::retain_fit_resources(resources, "caller_data"), "not owned")
  model_lease = dbreg:::retain_fit_resources(resources, names_owned[[1L]])
  output_lease = dbreg:::retain_fit_resources(resources, names_owned[[2L]])
  model = structure(list(resource_state = model_lease), class = "dbreg")
  dbreg:::release_fit_resources(resources)
  expect_equal(sort(DBI::dbListTables(conn)), sort(c("caller_data", "caller_export", names_owned[1:2])))
  expect_equal(DBI::dbGetQuery(dbreg:::fitted_resource_connection(model_lease), paste("SELECT * FROM", first))[[1L]], c(2, 4, 6))

  close(model)
  expect_error(dbreg:::fitted_resource_connection(model_lease), "released")
  expect_true(DBI::dbIsValid(conn))
  expect_equal(DBI::dbGetQuery(dbreg:::fitted_resource_connection(output_lease), paste("SELECT * FROM", second))[[1L]], c(2, 4, 6))
  # A repeated close must not drop a new caller table that reuses an old name.
  DBI::dbWriteTable(conn, names_owned[[1L]], data.frame(caller = 7))
  close(model)
  expect_equal(DBI::dbReadTable(conn, names_owned[[1L]])[[1L]], 7)
  dbreg:::fitted_resource_release(output_lease)
  expect_true(DBI::dbIsValid(conn))
  expect_true(DBI::dbExistsTable(conn, "caller_export"))

  # Serialized handles cannot be used as a live session or affect the original.
  next_table = manager[["create"]]("SELECT 9 AS learned")
  next_name = tail(resources[["tables"]], 1L)
  # An exhausted session cannot acquire new ownership; start a new fit manager.
  resources[["retained_session"]] = NULL
  lease = dbreg:::retain_fit_resources(resources, next_name)
  restored = unserialize(serialize(lease, NULL))
  expect_error(dbreg:::fitted_resource_connection(restored), "live owning session")
  dbreg:::fitted_resource_release(restored)
  expect_true(DBI::dbExistsTable(conn, next_name))
  expect_true(DBI::dbIsValid(dbreg:::fitted_resource_connection(lease)))
  dbreg:::fitted_resource_release(lease)

  # Finalization removes an abandoned owned table without closing a borrowed DB.
  resources[["retained_session"]] = NULL
  abandoned = manager[["create"]]("SELECT 5 AS learned")
  abandoned_name = tail(resources[["tables"]], 1L)
  lease = dbreg:::retain_fit_resources(resources, abandoned_name)
  rm(lease)
  invisible(gc())
  expect_false(DBI::dbExistsTable(conn, abandoned_name))
  expect_true(DBI::dbIsValid(conn))
  expect_error(close(model, unsupported = TRUE), "unsupported")

  # Portable models have nothing to release and keep predicting after close.
  portable = dbreg(y ~ x, data = data.frame(x = 1:8, y = c(2, 3, 5, 4, 6, 8, 7, 9)))
  before = predict(portable, data.frame(x = c(1, 4)))
  close(portable)
  expect_equal(predict(portable, data.frame(x = c(1, 4))), before)
})

local({
  # The final output lease keeps an owned connection alive after model close.
  resources = dbreg:::fit_resources()
  setup = dbreg:::setup_db_connection(NULL, NULL, data.frame(x = 1:3), NULL,
    resources = resources)
  conn = setup[["conn"]]
  on.exit(if (DBI::dbIsValid(conn)) DBI::dbDisconnect(conn, shutdown = TRUE))
  manager = dbreg:::fit_table_manager(conn, resources)
  first = manager[["create"]](paste("SELECT x*2 AS learned", setup[["from_statement"]]))
  second = manager[["create"]](paste("SELECT * FROM", first))
  model_lease = dbreg:::retain_fit_resources(resources, resources[["tables"]][[1L]])
  output_lease = dbreg:::retain_fit_resources(resources, resources[["tables"]][[1L]])
  dbreg:::release_fit_resources(resources)
  expect_true(DBI::dbIsValid(conn))
  expect_false(DBI::dbExistsTable(conn, setup[["registered_table"]]))
  dbreg:::fitted_resource_release(model_lease)
  expect_true(DBI::dbIsValid(conn))
  expect_equal(DBI::dbGetQuery(dbreg:::fitted_resource_connection(output_lease), paste("SELECT * FROM", second))[[1L]], c(2, 4, 6))
  dbreg:::fitted_resource_release(output_lease)
  expect_false(DBI::dbIsValid(conn))
  dbreg:::fitted_resource_release(output_lease)
})

local({
  conn = DBI::dbConnect(duckdb::duckdb())
  session = dbreg:::fitted_resource_session(conn)
  lease = dbreg:::fitted_resource_lease(session)
  DBI::dbDisconnect(conn, shutdown = TRUE)
  expect_error(dbreg:::fitted_resource_connection(lease), "live owning session")
  dbreg:::fitted_resource_release(lease)
  expect_error(dbreg:::fitted_resource_connection(lease), "released")
})

local({
  # A real dependency error must leave the failed drop owned and retryable,
  # while preserving the dependent relation created by the caller.
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbExecute(conn, "CREATE TABLE owned_parent(id INTEGER PRIMARY KEY)")
  DBI::dbExecute(conn, "CREATE TABLE caller_child(id INTEGER REFERENCES owned_parent(id))")
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  resources[["tables"]] = "owned_parent"
  lease = dbreg:::retain_fit_resources(resources, "owned_parent")
  expect_error(dbreg:::fitted_resource_release(lease), "Could not release")
  expect_true(DBI::dbExistsTable(conn, "caller_child"))
  expect_true(DBI::dbIsValid(dbreg:::fitted_resource_connection(lease)))
  DBI::dbExecute(conn, "DROP TABLE caller_child")
  dbreg:::fitted_resource_release(lease)
  expect_false(DBI::dbExistsTable(conn, "owned_parent"))
})

local({
  # Transient cleanup also retries failed drops and forgets successful ones.
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  DBI::dbExecute(conn, "CREATE TABLE transient_parent(id INTEGER PRIMARY KEY)")
  DBI::dbExecute(conn, "CREATE TABLE caller_child(id INTEGER REFERENCES transient_parent(id))")
  resources[["tables"]] = "transient_parent"
  expect_error(dbreg:::release_fit_resources(resources, quiet = FALSE), "Could not release fit resources")
  expect_identical(resources[["tables"]], "transient_parent")
  DBI::dbExecute(conn, "DROP TABLE caller_child")
  dbreg:::release_fit_resources(resources, quiet = FALSE)
  expect_identical(resources[["tables"]], character())
  DBI::dbWriteTable(conn, "transient_parent", data.frame(caller = 7))
  dbreg:::release_fit_resources(resources, quiet = FALSE)
  expect_equal(DBI::dbReadTable(conn, "transient_parent")[["caller"]], 7)
})
