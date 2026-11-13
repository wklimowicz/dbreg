# Native drivers may collect COUNT/BIGINT diagnostics as exact decimal text.
# Exercise real projection SQL while changing only the scalar transport mode.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  data = data.frame(g = c(1, 2, 3, 1, 2, 3), h = c(1, 2, 3, 2, 3, 1),
    w = c(1, 1, 1, 1e-12, 1e-12, 1e-12), v1 = c(1, -.2, 2, 5, 7, -1))
  DBI::dbWriteTable(conn, "qr_source", data)
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  transports = c("character", if (requireNamespace("bit64", quietly = TRUE)) "integer64")
  for (transport in transports) {
    project = dbreg:::dbglm_fe_qr_projection
    shim = new.env(parent = environment(project))
    shim$fe_cg_fetch = function(conn, sql, width, budget, limits) {
      result = dbreg:::fe_cg_fetch(conn, sql, width, budget, limits)
      result[] = lapply(result, if (transport == "character") as.character else bit64::as.integer64)
      result
    }
    environment(project) = shim
    run = function(limits) {
      budget = list2env(list(rows = 0, cells = 0, upload_bytes = 0), parent = emptyenv())
      project(conn, "qr_source", "v1", c("g", "h"),
        list(list(f = "g", v = "1.0"), list(f = "h", v = "1.0")),
        manager[["create"]], manager[["drop"]], budget, limits)
    }
    result = run(dbreg:::projection_limits())
    rows = DBI::dbGetQuery(conn, paste("SELECT * FROM", result[["table"]], "ORDER BY g,h"))
    X = cbind(outer(data$g, 1:3, `==`)*1, outer(data$h, 2:3, `==`)*1)
    reference = qr.resid(qr(X*sqrt(data$w)), data$v1*sqrt(data$w))/sqrt(data$w)
    expect_equal(rows$v1, reference[order(data$g, data$h)], tolerance = 1e-8)
    expect_identical(result[["method"]], "bounded_qr")
    dbreg:::release_fit_resources(resources, quiet = FALSE)
    expect_identical(DBI::dbListTables(conn), "qr_source")
    limited = dbreg:::projection_limits()
    limited[["memory_bytes"]] = 1
    expect_error(run(limited), "bounded QR recovery budget")
    expect_identical(DBI::dbListTables(conn), "qr_source")
  }
})
