# Controlled native-staging selection on DuckDB executes real CREATE/INSERT/
# DROP statements. Actual native JSON transport has a separate Databricks gate.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  resources = dbreg:::fit_resources()
  resources$conn = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  borrowed = manager$create("SELECT 17 AS keep")
  initial = sort(DBI::dbListTables(conn))
  owned = resources$tables
  set.seed(1493)
  pieces = list(matrix(rnorm(600), 200, 3), matrix(rnorm(450), 150, 3))
  lookup = dbreg:::fe_coefficient_lookups
  e = new.env(parent = environment(lookup))
  e$detect_backend = function(...) list(is_databricks = TRUE)
  environment(lookup) = e
  limits = dbreg:::projection_limits()
  budget = function() list2env(list(upload_bytes = 0), parent = emptyenv())
  actual = lookup(conn, pieces, budget(), limits, manager, statement_bytes = 4096)
  expect_equal(length(actual$tables), 2L)
  for (i in seq_along(pieces)) {
    rows = DBI::dbGetQuery(conn, paste("SELECT * FROM", actual$sources[i], "ORDER BY id"))
    expect_equal(rows$id, seq_len(nrow(pieces[[i]])))
    expect_equal(unname(as.matrix(rows[-1L])), pieces[[i]], tolerance = 1e-12)
  }
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 17)
  for (table in actual$tables) manager$drop(table)
  expect_identical(resources$tables, owned)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # Chunk overhead counts cumulatively, even when each family's full literal
  # would individually fit. A failed preflight must create no table.
  preview = budget()
  for (piece in pieces) invisible(dbreg:::fe_coefficient_sql(conn, piece, preview, limits))
  tight = limits
  tight$upload_bytes = preview$upload_bytes
  expect_error(lookup(conn, pieces, budget(), tight, manager, statement_bytes = 4096), "cumulative")
  expect_identical(resources$tables, owned)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  # A real INSERT schema failure occurs after a physical table was created.
  bad_manager = manager
  bad_manager$create = function(sql) manager$create(paste0("SELECT id,c1 FROM (", sql, ") incomplete"))
  expect_error(lookup(conn, pieces, budget(), limits, bad_manager, statement_bytes = 4096), "columns|values")
  expect_identical(resources$tables, owned)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 17)
})

local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  resources = dbreg:::fit_resources()
  resources$conn = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE, after = FALSE)
  manager = dbreg:::fe_cg_manager(conn, resources)
  borrowed = manager$create("SELECT 23 AS keep")
  set.seed(19175)
  d = expand.grid(g = 1:30, h = 1:20)
  d$x = rnorm(nrow(d))
  eta = .15 + .25 * d$x + .1 * sin(d$g) - .1 * cos(d$h)
  lookup = dbreg:::fe_coefficient_lookups
  e = new.env(parent = environment(lookup))
  e$detect_backend = function(...) list(is_databricks = TRUE)
  environment(lookup) = e
  step = dbreg:::dbglm_fe_step
  e = new.env(parent = environment(step))
  e$fe_coefficient_lookups = function(conn, pieces, budget, limits, manager) lookup(conn, pieces, budget, limits, manager, statement_bytes = 4096)
  environment(step) = e
  driver = dbreg:::dbglm_fe_irls
  e = new.env(parent = environment(driver))
  e$dbglm_fe_step = step
  environment(driver) = e
  control = dbglm_control(score_tol = 1e-9)
  for (family in list(binomial(), binomial("cloglog"))) {
    d$y = rbinom(nrow(d), 1L, family$linkinv(eta))
    base = data.frame(g = d$g, h = d$h, x = d$x)
    base[["__dbglm_y"]] = d$y
    base[["__dbglm_a"]] = 1
    base[["__dbglm_o"]] = 0
    DBI::dbWriteTable(conn, "base", base)
    fit = driver(conn, "base", "x", c("g", "h"), c(30L, 20L), nrow(d), 0,
      list(numeric(30), numeric(20)), list(link = family$link), control, lookup_manager = manager)
    dense = suppressWarnings(glm(y ~ x + factor(g) + factor(h), data = d, family = family,
      control = glm.control(epsilon = 1e-14, maxit = 1000L)))
    expected = unname(predict(dense))
    actual = as.numeric(d$x * fit$beta + fit$effects[[1L]][d$g] + fit$effects[[2L]][d$h])
    expect_true(fit$converged)
    expect_equal(as.numeric(fit$beta), coef(dense)[["x"]], tolerance = 1e-7)
    expect_equal(actual, expected, tolerance = 1e-7)
    expect_equal(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", fit$projected))$n, nrow(d))
    expect_equal(length(resources$tables), 3L) # borrowed marker + final two lookups
    expect_true(resources$peak_tables <= 5L) # Fisher + temporary Newton lookups
    expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 23)
    for (table in setdiff(as.character(DBI::dbQuoteIdentifier(conn, resources$tables)), borrowed)) manager$drop(table)
    DBI::dbRemoveTable(conn, "base")
  }
  # Warm projection also releases staged lookups when its identified basis
  # changes; original offsets remain outside this projection target.
  guess = eta + sin(seq_len(nrow(d))) / 5
  base = data.frame(g = d$g, h = d$h, x = d$x, copy = 2 * d$x)
  base[["__dbglm_a"]] = 1
  base[["__dbglm_o"]] = .1
  base[["__dbglm_etastart"]] = guess + .1
  DBI::dbWriteTable(conn, "base", base)
  warm = dbreg:::dbglm_warm_absorbed
  e = new.env(parent = environment(warm))
  e$dbglm_fe_step = step
  environment(warm) = e
  fit = warm(conn, "base", c("x", "copy"), c("g", "h"), c(30L, 20L), nrow(d),
    control, NULL, NULL, lookup_manager = manager)
  dense = lm(guess ~ x + factor(g) + factor(h), data = d)
  actual = as.numeric(d$x * fit$beta + fit$effects[[1L]][d$g] + fit$effects[[2L]][d$h])
  expect_identical(fit$keep, 1L)
  expect_equal(as.numeric(fit$beta), coef(dense)[["x"]], tolerance = 1e-9)
  expect_equal(actual, unname(fitted(dense)), tolerance = 1e-9)
  expect_equal(length(resources$tables), 1L)
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 23)
  DBI::dbRemoveTable(conn, "base")
  # Failed final certification rolls back newly staged lookups, not the marker.
  DBI::dbWriteTable(conn, "working", data.frame(g = d$g, h = d$h, x = d$x,
    w = .25, score = d$y - .5))
  e = environment(step)
  e$fe_projection_certificate = function(...) list(max_abs = Inf, max_relative = Inf)
  initial = sort(DBI::dbListTables(conn))
  owned = resources$tables
  expect_error(step(conn, "working", "x", c("g", "h"), c(30L, 20L), nrow(d), 1e-10,
    lookup_manager = manager), "certificate")
  expect_identical(resources$tables, owned)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  expect_equal(DBI::dbGetQuery(conn, paste("SELECT keep FROM", borrowed))$keep, 23)
  DBI::dbRemoveTable(conn, "working")
})
