# Shared native/local checks for compact sufficient-statistic transport.
dbreg_test_projection_transport = function(conn) {
  set.seed(4319)
  n = 1200L
  d = data.frame(id = seq_len(n), g = sample.int(13, n, TRUE), h = sample.int(9, n, TRUE),
    j = sample.int(5, n, TRUE), k = sample.int(3, n, TRUE), sex = sample(0:1, n, TRUE), part = sample(0:1, n, TRUE), w = runif(n, .2, 2))
  for (i in seq_len(10)) d[[paste0("x", i)]] = rnorm(n)
  d[["y"]] = .4*d[["x1"]] -.2*d[["x2"]] + sin(d[["g"]]) + rnorm(n)
  before = dbglm_validation_tables(conn)
  table = dbglm_validation_upload(conn, d)
  on.exit(dbglm_validation_drop(conn, table))
  fe = c("g", "h", "j", "k", "sex", "part")
  columns = c("y", paste0("x", seq_len(10)))
  terms = dbreg:::fe_value_terms(conn, columns)
  source = paste0("(SELECT * FROM ", table, ")")
  wide = dbreg:::read_fe_moments(conn, dbreg:::fe_moments_sql(conn, source, fe, terms))
  narrow = dbreg:::read_fe_moments(conn, dbreg:::fe_moments_sql(conn, source, fe, terms, compact = TRUE))
  a = dbreg:::assemble_fe_moments(dbreg:::decode_fe_moments(wide, length(fe), length(columns), 0L))
  b = dbreg:::assemble_fe_moments(dbreg:::decode_fe_moments(narrow, length(fe), length(columns), 0L))
  expect_equal(b[["B"]], a[["B"]], tolerance = 1e-13)
  expect_equal(b[["sys"]][["Cm"]], a[["sys"]][["Cm"]], tolerance = 1e-13)
  expect_equal(ncol(narrow), 4L)
  expect_true(is.null(narrow[["__b1"]]))
  sizes = vapply(d[fe], function(x) length(unique(x)), integer(1))
  required = 5*(nrow(narrow) + 1L) + sum(sizes)*length(columns)
  padded = nrow(wide)*ncol(wide)
  expect_true(required < padded)
  for (budget in c(required-1, required, required+1)) {
    summary = dbreg:::query_fe_moments(conn, source, fe, columns, max_cells = budget, fe_sizes = sizes)
    if (budget < required) expect_true(is.null(summary)) else {
      expect_equal(summary[["B"]], b[["B"]], tolerance = 1e-13)
      expect_true(summary[["transfer_cells"]] + 5 <= budget)
    }
  }
  # Fused GLM diagnostics have a sixth wire field and one global row. Their
  # budget must cover both payloads and the overflow probe, too.
  working = paste0("(SELECT *,0.0 AS z,0.0 AS eta,0.5 AS mu,1.0 AS dev,y AS score FROM ", table, ") working")
  fused_required = 6*(nrow(narrow)+2L)+sum(sizes)*3L+9L
  scalar = dbreg:::dbglm_fe_evaluate(conn, working, "x1", fe)
  for (budget in c(fused_required-1, fused_required, fused_required+1)) {
    limits = dbreg:::projection_limits()
    limits[["transfer_cells"]] = budget
    evaluated = dbreg:::dbglm_fe_evaluate(conn, working, "x1", fe,
      collect_fe_moments = TRUE, fe_sizes = sizes, limits = limits)
    expect_identical(is.null(evaluated[["fe_moments"]]), budget < fused_required)
    evaluated[["fe_moments"]] = NULL
    expect_equal(evaluated, scalar, tolerance = 1e-12)
  }
  old = options(dbreg.ap_solver = "moments", dbreg.ap_moment_max_cells = required,
    dbreg.ap_method = NULL, dbreg.duckdb_ap_method = NULL)
  on.exit(options(old), add = TRUE)
  formula = stats::as.formula(paste("y ~", paste(paste0("x", seq_len(10)), collapse = "+"), "|", paste(fe, collapse = "+")))
  fit = dbreg(formula, conn = conn, table = table, strategy = "demean", weights = "w", vcov = ~g)
  reference = fixest::feols(formula, d, weights = ~w, vcov = ~g,
    fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
  expect_identical(fit[["ap_method"]], "bounded_moments")
  expect_true(all(abs(coef(fit)-coef(reference)) <= 1e-8 + 1e-7*abs(coef(reference))))
  expect_true(all(abs(vcov(fit)-vcov(reference)) <= 1e-9 + 1e-6*abs(vcov(reference))))
  expect_identical(fit[["projection_diagnostics"]][[1L]][["status"]], "completed")
  expect_identical(dbglm_validation_tables(conn), sort(c(before, table)))
}
