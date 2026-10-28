# Forced-only discovery dictionaries must not outlive their local SQL source.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  index = seq_len(120)
  data = data.frame(id = index, f = rep(c("a", "b"), 60),
    x = sin(index * .19) + cos(index * .37), g = rep(seq_len(20), each = 6),
    w = 1 + (index %% 5) / 4)
  data[["y"]] = 1 + .6 * (data[["f"]] == "b") + 1.7 * data[["x"]] +
    sin(data[["g"]] * .23) + .03 * sin(index * .71)
  proof = Sys.getenv("DBREG_TEST_PROOF_DIR", "")
  for (case in c("forced", "ordinary_shared", "i_shared")) {
    learned = case == "forced"
    full_X = if (learned) function(rows) cbind(as.numeric(rows[["f"]] == "b"), rows[["x"]],
      outer(rows[["g"]], seq_len(20), `==`) * 1) else
      function(rows) cbind(1, rows[["x"]], as.numeric(rows[["f"]] == "b"))
    X = full_X(data)
    reference = stats::lm.wfit(X, data[["y"]], data[["w"]])
    expect_identical(reference[["rank"]], ncol(X))
    beta = reference[["coefficients"]]
    V = solve(crossprod(X, data[["w"]] * X)) *
      sum(data[["w"]] * reference[["residuals"]]^2) / (nrow(X) - ncol(X))
    formula = switch(case, forced = y ~ factor(f) + x | g,
      ordinary_shared = y ~ x + f + factor(f),
      i_shared = y ~ x + i(f, ref = "a") + factor(f))
    model = dbreg(formula, data = data, conn = conn, weights = "w",
      strategy = if (learned) "demean" else "moments",
      retain = if (learned) "prediction" else "none", ssc = list(K.fixef = "full"))
    on.exit(close(model), add = TRUE, after = FALSE)
    selected = if (learned) 1:2 else 1:3
    receipt = list(case = case, reference_beta = beta[selected], reference_full_V = V[selected, selected, drop = FALSE],
      actual_beta = coef(model), actual_full_V = vcov(model), stored_levels = model[["factor_levels"]])
    if (nzchar(proof)) { dir.create(proof, recursive = TRUE, showWarnings = FALSE); saveRDS(receipt, file.path(proof, paste0(case, "-before-prediction.rds"))) }
    expect_equal(as.numeric(receipt[["actual_beta"]]), as.numeric(receipt[["reference_beta"]]), tolerance = 1e-7)
    expect_equal(as.numeric(receipt[["actual_full_V"]]), as.numeric(receipt[["reference_full_V"]]), tolerance = 1e-7)
    rows = data[c(17, 2, 55, 9, 37, 79), ]
    rows[["f"]] = ifelse(rows[["f"]] == "a", "b", "a")
    rows[["x"]] = rows[["x"]] + .2
    rows[["y"]] = -999
    expected = as.numeric(full_X(rows) %*% beta)
    r = predict(model, rows)
    out = predict(model, rows, conn = conn, row_id = "id")
    table = as.data.frame(out)
    sql = table[[out[["prediction_column"]]]][match(rows[["id"]], table[["id"]])]
    close(out)
    design = dbreg:::sql_fitted_prediction_design(model, conn)
    receipt[["expected"]] = expected
    receipt[["R"]] = r
    receipt[["SQL"]] = sql
    receipt[["local_levels"]] = design[["factor_levels"]]
    receipt[["local_classes"]] = design[["factor_classes"]]
    receipt[["local_sources"]] = design[["sources"]]
    if (nzchar(proof)) saveRDS(receipt, file.path(proof, paste0(case, "-complete.rds")))
    expect_identical(is.na(r), is.na(expected))
    expect_equal(as.numeric(r), expected, tolerance = 1e-7)
    expect_identical(is.na(sql), is.na(expected))
    expect_equal(as.numeric(sql), expected, tolerance = 1e-7)
    expect_true(all(names(design[["factor_levels"]]) %in% names(design[["sources"]])))
    if (learned) {
      expect_false("f" %in% names(design[["factor_levels"]]))
      expect_false("f" %in% names(design[["factor_classes"]]))
      expect_true(all(c("f", "factor(f)") %in% names(model[["factor_levels"]])))
    } else expect_true("f" %in% names(design[["factor_levels"]]) && "f" %in% names(design[["sources"]]))
    bad = rows
    bad[["f"]][1L] = "unseen"
    expect_error(predict(model, bad), "unseen")
    close(model)
  }
})
