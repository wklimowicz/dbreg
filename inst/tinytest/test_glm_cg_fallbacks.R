# Large FE systems cannot hold the dense QR design. A projection that narrowly
# misses its certificate is then refined by CG, and weak working weights use
# certified CG instead of QR. Force both with a tiny design-entry budget: a
# loosened first CG target misses every certificate, and a quasi-separated
# zero spell (small prior weight) makes the working weights weak.
local({
  if (!requireNamespace("fixest", quietly = TRUE)) return(invisible())
  set.seed(11)
  nw = 30L; nt = 8L
  d = data.frame(g = rep(seq_len(nw), each = nt), t = rep(seq_len(nt), nw))
  home = sample(1:8, nw, TRUE)
  away = sample(1:8, nw, TRUE)
  d[["h"]] = ifelse(d[["t"]] > 4 & d[["g"]] %% 3 == 0, away[d[["g"]]], home[d[["g"]]])
  d[["x1"]] = rnorm(nrow(d))
  d[["x2"]] = rnorm(nrow(d))
  d[["y"]] = rpois(nrow(d), exp(-0.3 + 0.4 * d[["x1"]] - 0.2 * d[["x2"]] + 0.5 * rnorm(nw)[d[["g"]]]))
  d[["w"]] = 1
  # Workers 31 and 32 have positive outcomes at firm 9 only; worker 31's zero
  # spell at firm 3 is the only link to the main component (quasi-separation).
  b = data.frame(g = rep(31:32, each = nt), t = rep(seq_len(nt), 2L))
  spell = b[["g"]] == 31L & b[["t"]] > 4
  b[["h"]] = ifelse(spell, 3, 9)
  b[["x1"]] = rnorm(nrow(b))
  b[["x2"]] = rnorm(nrow(b))
  b[["y"]] = ifelse(spell, 0, rpois(nrow(b), 3) + 1)
  b[["w"]] = ifelse(spell, 1e-5, 1)
  d = rbind(d, b[names(d)])

  counts = new.env()
  counts[["qr"]] = counts[["declined"]] = counts[["refine"]] = 0
  qr = dbreg:::dbglm_fe_qr_projection
  cg = dbreg:::fe_sql_cg
  assignInNamespace("dbglm_fe_qr_projection", function(...) {
    counts[["qr"]] = counts[["qr"]] + 1
    tryCatch(qr(...), dbreg_projection_capacity = function(e) {
      counts[["declined"]] = counts[["declined"]] + 1
      stop(e)
    })
  }, "dbreg")
  # Only the step's first CG receives a warm-start argument; refinements do not.
  assignInNamespace("fe_sql_cg", function(conn, table, columns, identifiers, ops, scales, target, ..., start = NULL) {
    if (missing(start)) counts[["refine"]] = counts[["refine"]] + 1 else target = 1e4 * target
    cg(conn, table, columns, identifiers, ops, scales, target, ..., start = start)
  }, "dbreg")
  previous = options(dbreg.projection_max_entries = 400)
  on.exit({
    options(previous)
    assignInNamespace("dbglm_fe_qr_projection", qr, "dbreg")
    assignInNamespace("fe_sql_cg", cg, "dbreg")
  }, add = TRUE)

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  initial = sort(DBI::dbListTables(conn))
  fit = dbglm(y ~ x1 + x2 | g + h, conn = conn, data = d, weights = "w", family = poisson(),
    control = dbglm_control(projection_method = "sql"))
  expect_true(fit[["converged"]])
  expect_identical(fit[["projection_method"]], "sql_cg")
  # Every QR attempt was declined; some came from refinements, others from weak steps.
  expect_true(counts[["refine"]] > 0)
  expect_true(counts[["qr"]] > counts[["refine"]])
  expect_identical(counts[["declined"]], counts[["qr"]])

  reference = fixest::fepois(y ~ x1 + x2 | g + h, d, weights = ~w, notes = FALSE)
  expect_equal(nobs(fit), nobs(reference))
  expect_equal(coef(fit)[names(coef(reference))], coef(reference), tolerance = 1e-6)
  # The structural limit MLE drops the separated spell.
  limit = fixest::fepois(y ~ x1 + x2 | g + h, d[!(d[["g"]] == 31L & d[["t"]] > 4), ], weights = ~w, notes = FALSE)
  expect_equal(coef(fit)[names(coef(limit))], coef(limit), tolerance = 1e-8)
  close(fit)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
