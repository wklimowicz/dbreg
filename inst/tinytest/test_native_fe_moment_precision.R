# Evaluate the native aggregation expressions on a portable SQL engine. Packet
# transport is omitted here; actual Databricks receives its own integration gate.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  native = list(name = "spark", is_databricks = TRUE)
  summaries = dbreg:::fe_summary_sql
  e = new.env(parent = environment(summaries))
  e$detect_backend = function(...) native
  e$fe_packet_sql = function(conn, sql, ...) sql
  environment(summaries) = e
  certificate = dbreg:::fe_projection_certificate
  ce = new.env(parent = environment(certificate))
  ce$detect_backend = function(...) native
  environment(certificate) = ce
  DBI::dbExecute(conn, paste0("CREATE TEMP TABLE source AS SELECT i,1.0 AS w,",
    "CAST(i%2+1 AS INTEGER) AS g,CAST(i%3+1 AS INTEGER) AS k,",
    "CASE WHEN i<100000 THEN 1e5 WHEN i<200000 THEN 1e-7 ELSE -1e5 END AS v ",
    "FROM range(300000) r(i)"))
  # Independent integer counts determine the expected sums; no SQL SUM or
  # floating accumulation over the input is used as the reference.
  expected = list(c(.005, .005), c(1e5 + .0033333, .0033334, -1e5 + .0033333))
  for (ordering in c("i", "-i", "hash(i)")) {
    sql = summaries(conn, paste0("(SELECT * FROM source ORDER BY ", ordering, ") ordered"),
      c("g", "k"), "v", 20L)
    value = DBI::dbGetQuery(conn, sql)
    for (family in 1:2) {
      rows = value[value$kind == family, ]
      rows = rows[order(rows$i), ]
      expect_equal(rows$b1, expected[[family]], tolerance = 1e-10)
      expect_equal(rows$h, rep(300000 / nrow(rows), nrow(rows)))
    }
  }
  # This nonzero group mean must not disappear through cancellation and receive
  # a false orthogonality certificate. The two g groups have identical means.
  cert = certificate(conn, "source", "g", "v", 300000, 1e5)
  expect_equal(cert$max_abs, 1e-7 / 3, tolerance = 1e-12)
  expect_true(cert$max_abs > 1e-8)
  expect_true(is.finite(cert$max_relative))
  # Tiny values must retain their magnitude rather than being quantized to
  # zero by the guarded decimal branch. Large values take the finite fallback.
  for (units in c(1e-30, 1e12)) {
    DBI::dbExecute(conn, paste0("CREATE TEMP TABLE scaled AS SELECT g,k,w,",
      sprintf("%.17g", units), " AS v FROM source"))
    value = DBI::dbGetQuery(conn, summaries(conn, "scaled", c("g", "k"), "v", 20L))
    for (family in 1:2) {
      rows = value[value$kind == family, ]
      expect_true(all(is.finite(rows$b1) & rows$b1 != 0))
      expect_equal(rows$b1 / units, rep(300000 / nrow(rows), nrow(rows)), tolerance = 1e-10)
    }
    DBI::dbRemoveTable(conn, "scaled")
  }
  DBI::dbExecute(conn, "CREATE TEMP TABLE tiny_weights AS SELECT g,k,1e-30 AS w,1.0 AS v FROM source")
  tiny = DBI::dbGetQuery(conn, summaries(conn, "tiny_weights", c("g", "k"), "v", 20L))
  for (family in 1:2) {
    rows = tiny[tiny$kind == family, ]
    expect_true(all(is.finite(rows$h) & rows$h > 0))
    expect_equal(rows$h / 1e-30, rep(300000 / nrow(rows), nrow(rows)), tolerance = 1e-10)
  }
  DBI::dbRemoveTable(conn, "tiny_weights")
  # Decimal rounding must not hide a failing mean in a small-mass group.
  # The exact mean is 1e-10, above the unchanged 1e-11 certificate target.
  DBI::dbExecute(conn, "CREATE TEMP TABLE cancellation AS SELECT 1 AS g,2e-9 AS w,1.0 AS v UNION ALL SELECT 1,2e-9,-1.0+2e-10")
  small = certificate(conn, "cancellation", "g", "v", 2L, 1)
  expect_true(small$max_abs > 1e-11)
  expect_true(small$max_relative > 1e-11)
  DBI::dbRemoveTable(conn, "cancellation")
  for (mass in c(2e-9, 1e-30, 1)) {
    DBI::dbExecute(conn, paste0("CREATE TEMP TABLE balanced AS SELECT 1 AS g,",
      sprintf("%.17g", mass), " AS w,1.0 AS v UNION ALL SELECT 1,",
      sprintf("%.17g", mass), ",-1.0"))
    zero = certificate(conn, "balanced", "g", "v", 2L, 1)
    expect_equal(zero$max_abs, 0)
    expect_equal(zero$max_relative, 0)
    DBI::dbRemoveTable(conn, "balanced")
  }
  DBI::dbRemoveTable(conn, "source")
  expect_identical(DBI::dbListTables(conn), character())
})

# A classified failure at an already unchanged initial anchor cannot improve
# through coefficient halving. Preserve its condition rather than repeating SQL.
local({
  driver = dbreg:::dbglm_fe_irls
  e = new.env(parent = environment(driver))
  state = new.env(parent = emptyenv())
  state$attempts = 0L
  e$dbglm_fe_trial = function(...) list(bad = 0L, mu_min = .5, mu_max = .5,
    deviance = 10, H = matrix(1, 1, 1), score = 0, working = "mock")
  e$dbglm_fe_step = function(...) {
    state$attempts = state$attempts + 1L
    dbreg:::projection_abort("certificate_failed", "unchanged-anchor certificate",
      list(max_abs = 1e-8, max_relative = 1e-7, target = 1e-10))
  }
  environment(driver) = e
  result = tryCatch(driver(NULL, "mock", "x", "g", 1L, 10L, 0, list(0),
    list(link = "identity"), dbglm_control(step_max_halvings = 25L)), error = identity)
  expect_true(inherits(result, "dbreg_projection_numerical"))
  expect_identical(result$reason, "certificate_failed")
  expect_equal(result$details$target, 1e-10)
  expect_equal(state$attempts, 1L)
})
