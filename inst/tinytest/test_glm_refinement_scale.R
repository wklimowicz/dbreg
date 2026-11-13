# Weighted row replication changes moment magnitudes, not the within projection.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = expand.grid(g = 1:7, h = 1:5, rep = 1:10)
  i = seq_len(nrow(d))
  d$x = sin(i / 7) + .01 * d$g + .03 * d$h
  response = cos(i / 3) + .2 * d$g - .1 * d$h
  a = rep(c(.5, 1.5, 4), length.out = nrow(d))
  D = model.matrix(~ factor(g) + factor(h), d)
  Y = cbind(response, d$x)
  dense = qr.solve(D * sqrt(a), Y * sqrt(a))
  expected = Y - D %*% dense
  project = dbreg:::dbglm_fe_step
  for (mass in c(1, 1e20)) {
    working = transform(d, w = a * mass, score = a * mass * response)
    DBI::dbWriteTable(conn, "working", working)
    fit = tryCatch(project(conn, "working", "x", c("g", "h"), c(7L, 5L), nrow(d), 1e-10), error = identity)
    expect_false(inherits(fit, "error"))
    if (!inherits(fit, "error")) {
      actual = Y - fit$projection[[1L]][d$g, , drop = FALSE] - fit$projection[[2L]][d$h, , drop = FALSE]
      expect_equal(unname(actual), unname(expected), tolerance = 1e-8)
      sql = DBI::dbGetQuery(conn, paste("SELECT SUM(w*v1*v1) AS a,SUM(w*v1*v2) AS b,SUM(w*v2*v2) AS c FROM", fit$residual_source))
      gram = crossprod(expected, expected * a)
      expect_equal(as.numeric(sql[1, ]) / mass, c(gram[1, 1], gram[1, 2], gram[2, 2]), tolerance = 1e-7)
      expect_true(fit$certificate$max_abs <= 1e-10 && fit$certificate$max_relative <= 1e-10)
      expect_true(is.finite(fit$refinement$residual))
    }
    DBI::dbRemoveTable(conn, "working")
    expect_identical(DBI::dbListTables(conn), character())
  }
  DBI::dbWriteTable(conn, "working", transform(d, w = a, score = a * response))
  # Norm collection is charged before marginal/pair summaries are requested.
  e = new.env(parent = environment(project))
  captured = new.env(parent = emptyenv())
  captured$recorded = NULL
  e$query_fe_summaries = function(conn, source, fe, columns, sizes, budget, limits) {
    captured$recorded = list(cells = budget$cells, rows = budget$rows)
    stop("stop after norm dispatch")
  }
  environment(project) = e
  expect_error(project(conn, "working", "x", c("g", "h"), c(7L, 5L), nrow(d), 1e-10), "after norm dispatch")
  expect_equal(captured$recorded, list(cells = 6L, rows = 1L))
  prior = options(dbreg.projection_max_cells = 3L)
  expect_error(project(conn, "absent", "x", c("g", "h"), c(7L, 5L), nrow(d), 1e-10), "working-column norms exceed")
  options(prior)
  environment(project) = asNamespace("dbreg")
  # An inner roundoff-floor exit never bypasses the final SQL row certificate.
  e = new.env(parent = environment(project))
  captured$floor_seen = FALSE
  e$dbglm_fe_refine = function(...) {
    value = dbreg:::dbglm_fe_refine(...)
    captured$floor_seen = value$refinement$roundoff_floor
    value
  }
  e$fe_projection_certificate = function(...) list(max_abs = Inf, max_relative = Inf)
  environment(project) = e
  expect_error(project(conn, "working", "x", c("g", "h"), c(7L, 5L), nrow(d), 1e-25), "certificate")
  expect_true(captured$floor_seen)
  expect_identical(DBI::dbListTables(conn), "working")
  DBI::dbRemoveTable(conn, "working")
  DBI::dbWriteTable(conn, "working", transform(d, x = 0, w = a, score = 0))
  zero = dbreg:::dbglm_fe_step(conn, "working", "x", c("g", "h"), c(7L, 5L), nrow(d), 1e-10)
  expect_equal(as.numeric(zero$H), 0)
  expect_equal(as.numeric(zero$score), 0)
  expect_identical(zero$refinement$roundoff_floor, FALSE)
  DBI::dbRemoveTable(conn, "working")
  expect_identical(DBI::dbListTables(conn), character())
})

# Isolate the stopping policy with controlled equation roundoff; real SQL and
# dense-reference projection checks above remain separate correctness evidence.
local({
  probe = function(errors, mass = 1, units = 1, tol = 1e-20, max_iter = 8L) {
    refine = dbreg:::dbglm_fe_refine
    e = new.env(parent = environment(refine))
    state = new.env(parent = emptyenv())
    state$call = 0L
    e$fe_multiply = function(system, projection) {
      state$call = state$call + 1L
      value = errors[min(state$call, length(errors))]
      matrix(-value * sqrt(mass) * abs(units) * sqrt(mass), 1, 1)
    }
    environment(refine) = e
    solver = list(solve = function(rhs) matrix(0, 1, 1))
    refine(list(diag = mass), matrix(0, 1, 1), solver,
      abs(units) * sqrt(mass), tol, max_iter)
  }
  for (mass in c(1e-10, 1, 1e10)) for (units in c(1e-10, 1, 1e10)) {
    fit = probe(c(1e-15, 1e-15), mass, units)
    expect_true(fit$refinement$roundoff_floor)
    expect_equal(fit$refinement$residual / 1e-15, 1, tolerance = 1e-10)
  }
  expect_false(probe(1e-15, tol = 1e-10)$refinement$roundoff_floor)
  expect_error(probe(1e-6), "did not converge")
  expect_error(probe(c(4e-14, 1e-14, 2e-15), max_iter = 3L), "did not converge")
  expect_error(probe(1e-15, max_iter = 1L), "did not converge")
  expect_error(probe(Inf), "did not converge")
  expect_error(probe(NA_real_), "did not converge")
})
