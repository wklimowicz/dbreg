local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  set.seed(11146)
  d = expand.grid(g = 1:5, h = 1:4, k = 1:3)
  d[["w"]] = runif(nrow(d), .25, 2)
  d[["x"]] = rnorm(nrow(d)) + d[["g"]]
  d[["z"]] = rnorm(nrow(d)) - d[["k"]]
  d[["row_id"]] = seq_len(nrow(d))
  DBI::dbWriteTable(conn, "cg_fit", d)
  fe = c("g", "h", "k")
  X = model.matrix(~ factor(g) + factor(h) + factor(k), d)
  oracle = lapply(d[c("x", "z")], function(v) lm.wfit(X, v, d[["w"]])[["residuals"]])
  fit = dbreg:::fe_sql_cg(conn, "cg_fit", c("x", "z"), c("row_id", fe),
    lapply(fe, function(f) list(f = f, v = "1.0")), c(1, 1), 1e-12,
    dbreg:::fe_cg_control(), create, drop,
    list2env(list(rows = 0, cells = 0), parent = emptyenv()),
    dbreg:::projection_limits(), retain_coefficients = TRUE)
  coefs = DBI::dbGetQuery(conn, paste("SELECT * FROM", fit[["coefficients"]]))
  residuals = DBI::dbGetQuery(conn, paste("SELECT x,z FROM", fit[["table"]], "ORDER BY row_id"))
  expect_equal(nrow(coefs), sum(vapply(d[fe], function(v) length(unique(v)), integer(1))))
  for (i in seq_along(oracle)) {
    adjustment = Reduce(`+`, lapply(seq_along(fe), function(j) {
      family = coefs[coefs[["__cg_op"]] == j, ]
      family[[paste0("c", i)]][match(d[[fe[j]]], family[["__cg_group"]])]
    }))
    expect_equal(d[[names(oracle)[i]]] - adjustment, oracle[[i]], tolerance = 1e-10)
    expect_equal(residuals[[i]], oracle[[i]], tolerance = 1e-10)
  }
  drop(fit[["coefficients"]])
  drop(fit[["table"]])
  expect_identical(DBI::dbListTables(conn), "cg_fit")
})
