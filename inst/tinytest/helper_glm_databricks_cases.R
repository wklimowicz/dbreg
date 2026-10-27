# Run the same numerical cases through the local SQL fixture and a live
# Databricks connection. References use only fixest's own estimates/refinements.
dbglm_test_databricks_cases = function(conn, verbose = FALSE) {
  # Wide FE covariance queries embed large coefficient literals. Verify their
  # translation and execution on both the fixture and the native warehouse.
  payload = strrep("QUJDZGVmMDEyMzQ1Njc4OSsv", 300000L)
  decoded = dbreg:::backend_get_query(conn, paste0("SELECT LENGTH('", payload,
    "') AS n,", dbreg:::sql_sum(conn, "x"), " AS total FROM (SELECT 1.0 AS x UNION ALL SELECT 2.0 AS x) t"))
  expect_equal(as.numeric(decoded[["n"]]), nchar(payload))
  expect_equal(as.numeric(decoded[["total"]]), 3)
  # A wide coefficient dictionary must be shared across FE operations. Six
  # copies of this bounded lookup exceeded the real warehouse query limit.
  coefficients = matrix(sin(seq_len(2400L * 88L)), ncol = 88L)
  info = rep(list(list(id = seq_len(400L))), 6L)
  frame = dbreg:::dbglm_coefficient_frame(info, coefficients)
  source = paste0("(", dbreg:::sql_values_frame(conn, frame), ")")
  ctes = dbreg:::dbglm_coefficient_ctes(source, 6L, "wide")
  branches = vapply(seq_len(6L), function(i) paste0("SELECT ", i,
    " AS op,COUNT(*) AS n,SUM(c1) AS first,SUM(c88) AS last FROM wide_g", i), character(1))
  query = paste0("WITH ", paste(ctes, collapse = ","), " ", paste(branches, collapse = " UNION ALL "), " ORDER BY op")
  expect_true(nchar(query, type = "bytes") < 8 * 1024^2)
  stored = dbreg:::backend_get_query(conn, query)
  expect_identical(as.numeric(stored[["op"]]), as.numeric(seq_len(6L)))
  expect_equal(as.numeric(stored[["n"]]), rep(400, 6L))
  expect_equal(stored[["first"]], vapply(split(coefficients[, 1], frame[["__glm_op"]]), sum, numeric(1)),
    tolerance = 1e-11, check.attributes = FALSE)
  expect_equal(stored[["last"]], vapply(split(coefficients[, 88], frame[["__glm_op"]]), sum, numeric(1)),
    tolerance = 1e-11, check.attributes = FALSE)
  owned = new.env(parent = emptyenv())
  owned[["tables"]] = character()
  on.exit(for (table in owned[["tables"]]) try(DBI::dbExecute(conn, paste("DROP TEMP TABLE", table)), silent = TRUE))
  upload = function(data) {
    table = basename(tempfile("glm_test_"))
    owned[["tables"]] = c(owned[["tables"]], table)
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", table, "AS", dbreg:::sql_values_frame(conn, data)))
    table
  }
  compare = function(fit, reference) {
    expect_identical(names(coef(fit)), names(coef(reference)))
    expect_true(all(abs(coef(fit) - coef(reference)) <= 1e-7 + 1e-6 * abs(coef(reference))))
    expect_true(all(abs(vcov(fit) - vcov(reference)) <= 1e-8 + 1e-5 * abs(vcov(reference))))
    expect_true(all(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(reference)))) <= 1e-8 + 1e-6 * sqrt(diag(vcov(reference)))))
    expect_equal(nobs(fit), nobs(reference))
    expect_true(abs(as.numeric(logLik(fit) - logLik(reference))) / nobs(fit) < 1e-9)
    expect_equal(fit[["n_params_model"]], reference[["nparams"]])
    expect_equal(fit[["n_params_vcov"]], attr(vcov(reference, attr = TRUE), "df.K"))
    expect_true(fit[["converged"]] && fit[["score_max"]] <= fit[["control"]][["score_tol"]])
    expect_true(all(diff(fit[["deviance_history"]]) <= 64 * .Machine$double.eps * pmax(1, head(fit[["deviance_history"]], -1))))
  }
  set.seed(7235)
  n = 1600L
  d = data.frame(id = seq_len(n), f1 = rep(seq_len(20), each = 80),
    f2 = sample(seq_len(9), n, TRUE), f3 = sample(seq_len(7), n, TRUE),
    f4 = sample(seq_len(5), n, TRUE), f5 = rbinom(n, 1, .7), f6 = rbinom(n, 1, .25),
    a = sample(1:3, n, TRUE), o = runif(n, -.15, .15), u = rnorm(n), z = rnorm(n))
  # Deliberate dependence between FE dimensions and within-FE variation in x.
  d[["f2"]][seq_len(1000)] = (d[["f1"]][seq_len(1000)] - 1L) %% 9L + 1L
  d[["x"]] = 10.8 + .4 * sin(d[["f1"]]) + .1 * d[["u"]]
  eta = -1 + .3 * d[["u"]] - .2 * d[["z"]] + .2 * sin(d[["f1"]]) + d[["o"]]
  formula = y ~ x + z | f1 + f2 + f3 + f4 + f5 + f6
  for (family in list(binomial("cloglog"), binomial("logit"), binomial("probit"), poisson())) {
    d[["y"]] = if (family[["family"]] == "poisson") rpois(n, exp(eta)) else rbinom(n, 1, family[["linkinv"]](eta))
    table = upload(d)
    for (kind in c("iid", "hc1", "cluster")) {
      if (verbose) message("[Databricks test] ", family[["family"]], "/", family[["link"]], " ", kind)
      vc = if (kind == "cluster") ~f2 else kind
      ref = dbglm_reference_feglm(formula, d, family = family, weights = ~a, offset = ~o,
        vcov = if (kind == "hc1") "hetero" else vc, fixef.rm = "none")
      fit = dbglm(formula, conn = conn, table = table, family = family, weights = "a", offset = "o", vcov = vc)
      compare(fit, ref)
      expect_identical(fit[["projection_method"]], "moments")
      expect_equal(unname(fit[["design_center"]]), colMeans(d[c("x", "z")]), check.attributes = FALSE, tolerance = 1e-12)
    }
    if (family[["link"]] == "cloglog") {
      # The join fallback must agree with arrays under the same estimator.
      saved = options(dbreg.glm_databricks_array_max_entries = 0)
      fallback = tryCatch(dbglm(formula, conn = conn, table = table, family = family,
        weights = "a", offset = "o", vcov = ~f2), finally = options(saved))
      compare(fallback, ref)
      expect_equal(coef(fallback), coef(fit), tolerance = 1e-8)
      expect_equal(vcov(fallback), vcov(fit), tolerance = 1e-8)
      # Translating an entire design column by a large constant must not
      # change which slopes are identified or the full covariance matrix.
      shifted = d
      shifted[["x"]] = shifted[["x"]] + 1e6
      shifted_table = upload(shifted)
      shifted_fit = dbglm(formula, conn = conn, table = shifted_table, family = family,
        weights = "a", offset = "o", vcov = ~f2)
      compare(shifted_fit, ref)
      expect_identical(names(coef(shifted_fit)), c("x", "z"))
      expect_equal(coef(shifted_fit), coef(fit), tolerance = 1e-7)
      expect_equal(vcov(shifted_fit), vcov(fit), tolerance = 1e-7)
    }
    expect_equal(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]]), n)
    expect_error(dbglm(formula, conn = conn, table = table, family = family,
      control = dbglm_control(glm_max_iter = 1)), "converge")
    expect_true(DBI::dbIsValid(conn))
  }

  # Check missing/sparse rank semantics with actual numerical results, not SQL
  # substrings. This also exercises double decoding on a real warehouse.
  ranks = data.frame(f1 = c(1, 2, 3), f2 = c(1, 2, 1), v = c(2, -1, .5))
  table = upload(ranks)
  ops = list(list(f = "f1", v = "1.0"), list(f = "f2", v = "v"))
  info = list(list(id = c(1, 3)), list(id = c(1, 2)))
  coefficients = c(.12345678901234567, -2.5, .7, -.1)
  terms = dbreg:::dbglm_databricks_array_terms(ops, info, coefficients, conn = conn)
  value = DBI::dbGetQuery(conn, paste0("SELECT ", paste(terms, collapse = " + "), " AS eta FROM ", table, " s ORDER BY f1"))[["eta"]]
  expect_equal(value, c(coefficients[[1]] + 1.4, .1, -2.5 + .35), tolerance = 1e-14)
  sparse = info
  sparse[[1]][["id"]] = c(1, 1000001)
  expect_null(dbreg:::dbglm_databricks_array_terms(ops, sparse, coefficients, conn = conn))
  expect_null(dbreg:::dbglm_databricks_array_terms(ops, info, coefficients, max_bytes = 1, conn = conn))
  expect_null(dbreg:::dbglm_databricks_array_terms(ops, info, coefficients, max_entries = 0, conn = conn))
  expect_error(dbreg:::dbglm_databricks_array_terms(ops, info, c(NA, coefficients[-1]), conn = conn), "Nonfinite")

  # A bounded fused query must fall back to complete scalar diagnostics when
  # even the FE marginals cannot fit, never certify a truncated result.
  working = data.frame(f1 = c(1, 1, 2, 2), f2 = c(1, 2, 1, 2), x = c(1, 2, 3, 4),
    z = c(2, 0, 1, -1), w = c(1, 2, 3, 4), eta = 0, mu = .5, dev = 1, score = c(1, -2, 3, -4))
  table = upload(working)
  fe = c("f1", "f2"); attr(fe, "slopes") = list(NULL, NULL)
  full = dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = FALSE)
  for (cap in c(1, 3, 7)) {
    saved = options(dbreg.duckdb_ap_max_edges = cap)
    capped = tryCatch(dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = TRUE), finally = options(saved))
    expect_null(capped[["fe_moments"]])
    expect_equal(capped, full, tolerance = 1e-13)
  }
  # Overflow can omit the global row entirely: reject it before JSON decoding.
  expect_null(dbreg:::read_fe_moments(conn,
    "SELECT 1 AS __set,CAST(NULL AS STRING) AS __globals", max_rows = 0))
  fused = dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = TRUE)
  expect_equal(fused[["fe_moments"]][["rows"]], 8)
  cached = dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = TRUE,
    moment_rows = fused[["fe_moments"]][["rows"]])
  expect_equal(cached, fused, tolerance = 1e-13)
  expect_error(dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = TRUE,
    moment_rows = 7), "support changed")
  expect_error(dbreg:::dbglm_fe_evaluate(conn, table, "x", fe, collect_fe_moments = TRUE,
    moment_rows = 9), "support changed")
  saved = options(dbreg.duckdb_ap_max_edges = 1)
  recapped = tryCatch(dbreg:::dbglm_fe_evaluate(conn, table, "x", fe,
    collect_fe_moments = TRUE, moment_rows = 8), finally = options(saved))
  expect_equal(recapped, full, tolerance = 1e-13)
  expect_equal(fused[["raw_moments"]], unname(crossprod(as.matrix(working[c("x", "z")]), working[["w"]] * as.matrix(working[c("x", "z")]))))
  expect_equal(fused[["score_max"]], full[["score_max"]], tolerance = 1e-13)
  expect_equal(fused[["deviance"]], 4)
  compact = dbreg:::read_fe_moments(conn, dbreg:::fe_moments_sql(conn,
    table, fe, dbreg:::fe_value_terms(conn, c("x", "z"), "score"),
    c("SUM(dev) AS deviance", "SUM(score) AS raw_score")))
  expect_equal(attr(compact, "dbglm_globals"), data.frame(deviance = 4, raw_score = -2))
  # Resource contract: a global statistic must occupy one scalar, rather than
  # being padded over all FE-pair rows as the regressor count grows.
  expect_false("deviance" %in% names(compact))
  extremes = dbreg:::read_fe_moments(conn, paste0("SELECT 0 AS __set,",
    "to_json(named_struct('large',a,'small',b,'fraction',c)) AS __globals FROM ",
    "(SELECT CAST(1e300 AS DOUBLE) AS a,CAST(1e-300 AS DOUBLE) AS b,",
    "CAST(0.12345678901234567 AS DOUBLE) AS c) t"))
  direct = DBI::dbGetQuery(conn, paste0("SELECT CAST(1e300 AS DOUBLE) AS a,",
    "CAST(1e-300 AS DOUBLE) AS b,CAST(0.12345678901234567 AS DOUBLE) AS c"))
  # Compare transport with the database's own doubles: decimal-to-double
  # conversion can differ from R's literal parser by one rounding unit.
  transported = as.numeric(attr(extremes, "dbglm_globals")[1, ])
  database_doubles = as.numeric(direct[1, ])
  if (!identical(transported, database_doubles)) {
    print(data.frame(field = names(attr(extremes, "dbglm_globals")),
      transported = sprintf("%a", transported),
      database = sprintf("%a", database_doubles),
      relative_error = (transported - database_doubles) / database_doubles))
  }
  expect_identical(as.numeric(attr(extremes, "dbglm_globals")[1, ]), as.numeric(direct[1, ]))
  expect_true(all(abs(as.numeric(direct[1, ]) / c(1e300, 1e-300, .12345678901234567) - 1) <= 4 * .Machine$double.eps))
  expect_error(dbreg:::read_fe_moments(conn,
    "SELECT 1 AS __set,CAST(NULL AS STRING) AS __globals"), "global summary")
  expect_error(dbreg:::validate_moment_limits(NA, 1), "limits")
  expect_error(dbreg:::validate_moment_limits(1, Inf), "limits")
  expect_error(dbreg:::validate_moment_limits(.5, 1), "limits")
}
