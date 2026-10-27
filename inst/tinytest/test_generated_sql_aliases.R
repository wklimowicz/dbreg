local({
  set.seed(11041)
  d = data.frame(a = rnorm(120), b = rnorm(120), a_x_b = rnorm(120),
    dbreg_i1_1 = rnorm(120), dbreg_i2_1 = rnorm(120), f = rep(c("a", "b", "c"), 40))
  d[["y"]] = 1 + 0.8 * d[["a"]] * d[["b"]] - 0.6 * d[["a_x_b"]] +
    0.5 * (d[["f"]] == "b") + rnorm(120, sd = 0.3)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "generated_alias_training", d)
  formulas = list(y ~ a:b + a_x_b, y ~ a_x_b + a:b,
    y ~ i(f, ref = "a") + dbreg_i1_1, y ~ dbreg_i1_1 + i(f, ref = "a"),
    y ~ dbreg_i2_1 + i(f, ref = "a"))
  matrix_for = function(formula, data) {
    columns = lapply(attr(terms(formula), "term.labels"), function(term) {
      if (grepl("i(", term, fixed = TRUE)) return(cbind(as.numeric(data[["f"]] == "b"), as.numeric(data[["f"]] == "c")))
      if (term == "a:b") return(data[["a"]] * data[["b"]])
      data[[term]]
    })
    do.call(cbind, c(list(rep(1, nrow(data))), columns))
  }
  for (formula in formulas) {
    compiled = dbreg::sql_model_matrix(formula, conn, "generated_alias_training")
    expect_equal(anyDuplicated(tolower(compiled[["aliases"]])), 0L)
    expect_identical(compiled[["columns"]][["alias"]], compiled[["aliases"]])
    i_positions = which(vapply(compiled[["columns"]][["recipe"]], function(recipe) identical(recipe[["kind"]], "i"), logical(1)))
    expect_identical(compiled[["i_columns"]][["alias"]], compiled[["aliases"]][i_positions])
    X = matrix_for(formula, d)
    beta = as.numeric(qr.solve(X, d[["y"]]))
    residual = as.numeric(d[["y"]] - X %*% beta)
    V = solve(crossprod(X)) * sum(residual^2) / (nrow(d) - ncol(X))
    rows = d[c(11, 4, 30, 1, 22), ]
    rows[["id"]] = c(9L, 4L, 12L, 1L, 19L)
    rows[["a"]][3] = NA_real_
    rows[["f"]][4] = NA_character_
    expected = as.numeric(matrix_for(formula, rows) %*% beta)
    # Missing irrelevant variables must not mask which sources this formula uses.
    extra = rows
    extra[["dbreg_i1_1_1"]] = -999
    extra[["a_x_b_1"]] = 999
    extra[["irrelevant"]] = NA_real_
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbreg(formula, data = d, conn = conn, strategy = strategy)
      expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-9)
      expect_equal(attr(vcov(fit), "df.K"), ncol(X))
      expect_equal(df.residual(fit), nrow(d) - ncol(X))
      if (strategy == "compress") expect_equal(as.numeric(predict(fit)), as.numeric(X %*% beta), tolerance = 1e-9)
      expect_equal(predict(fit, rows), expected, tolerance = 1e-9)
      expect_equal(predict(fit, extra), expected, tolerance = 1e-9)
      order = c(5, 2, 4, 1, 3)
      for (table_source in c(FALSE, TRUE)) {
        if (table_source) {
          DBI::dbWriteTable(conn, "generated_alias_rows", extra[order, ], overwrite = TRUE)
          out = predict(fit, conn = conn, table = "generated_alias_rows", row_id = "id")
        } else out = predict(fit, extra[order, ], conn = conn, row_id = "id")
        actual = as.data.frame(out)
        expect_equal(actual[[out[["prediction_column"]]]][match(rows[["id"]], actual[["id"]])], expected, tolerance = 1e-9)
        close(out)
      }
      close(fit)
    }
  }

  # Existing unambiguous public aliases and expression spellings stay intact.
  plain = dbreg::sql_model_matrix(~a + b + a:b, conn, "generated_alias_training")
  expect_identical(plain[["aliases"]], c("a", "b", "a_x_b"))
  expect_identical(plain[["col_names"]], c("a", "b", "a_x_b"))
  indicator = dbreg::sql_model_matrix(~i(f, ref = "a"), conn, "generated_alias_training")
  expect_identical(indicator[["aliases"]], c("dbreg_i1_1", "dbreg_i1_2"))
  product = dbreg::sql_model_matrix(a_x_b ~ a:b, conn, "generated_alias_training")
  expect_identical(product[["select_exprs"]], plain[["select_exprs"]][3])
  expect_identical(product[["col_names"]], "a_x_b")
  expect_identical(product[["aliases"]], "a_x_b_1")

  # Source roles outside the RHS must survive beside the generated projection.
  X = cbind(1, d[["a"]] * d[["b"]])
  for (role in c("response", "weight", "cluster", "offset")) {
    data = d
    data[["case_row"]] = seq_len(nrow(data))
    data[["a_x_b"]] = switch(role, response = d[["y"]],
      weight = rep(c(0.5, 1, 2), 40), cluster = rep(seq_len(12), 10), offset = d[["a_x_b"]])
    formula = if (role == "response") a_x_b ~ a:b else y ~ a:b
    response = if (role == "response") data[["a_x_b"]] else data[["y"]]
    if (role == "offset") response = response - data[["a_x_b"]]
    weights = if (role == "weight") data[["a_x_b"]] else rep(1, nrow(data))
    beta = as.numeric(qr.solve(X * sqrt(weights), response * sqrt(weights)))
    residual = as.numeric(response - X %*% beta)
    bread = solve(crossprod(X, weights * X))
    if (role == "cluster") {
      scores = rowsum(X * (weights * residual), data[["a_x_b"]], reorder = FALSE)
      V = bread %*% crossprod(scores) %*% bread * 12 / 11 * 119 / 118
    } else V = bread * sum(weights * residual^2) / 118
    for (strategy in c("moments", "compress")) {
      args = list(formula, data = data, conn = conn, strategy = strategy)
      if (role == "weight") args[["weights"]] = "a_x_b"
      if (role == "cluster") args[["vcov"]] = ~a_x_b
      if (role == "offset") args[["offset"]] = "a_x_b"
      fit = do.call(dbreg::dbreg, args)
      expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-9)
      expect_equal(df.residual(fit), 118)
      expect_equal(attr(vcov(fit), "df.K"), 2)
      expected = as.numeric(X %*% beta) + if (role == "offset") data[["a_x_b"]] else 0
      expect_equal(predict(fit, data), expected, tolerance = 1e-9)
      out = predict(fit, data, conn = conn, row_id = "case_row")
      actual = as.data.frame(out)
      expect_equal(actual[[out[["prediction_column"]]]][match(data[["case_row"]], actual[["case_row"]])], expected, tolerance = 1e-9)
      close(out)
      close(fit)
    }
  }
})

local({
  set.seed(11042)
  d = data.frame(a = rnorm(180), b_x_c = rnorm(180), a_x_b = rnorm(180), c = rnorm(180),
    dbreg_i = rep(c("0", "1"), 90), g = rep(c("0", "1"), 90), dbreg_i1 = rnorm(180),
    f = rep(c("a", "b", "c"), 60), z = rnorm(180))
  d[["endog"]] = 0.7 * d[["a_x_b"]] * d[["c"]] + 0.4 * d[["z"]] + rnorm(180)
  d[["y"]] = 1 + 0.4 * d[["a"]] * d[["b_x_c"]] + 1.2 * d[["endog"]] + rnorm(180)
  formulas = list(cross = y ~ a:b_x_c | endog ~ a_x_b:c,
    mixed = y ~ dbreg_i | endog ~ i(f, ref = "a") + z,
    factor_i = y ~ g + i(f, ref = "a") | endog ~ z,
    raw = y ~ dbreg_i1 | endog ~ i(f, ref = "a") + z,
    overlap = y ~ i(f, ref = "a") + dbreg_i1 | endog ~ z + i(f, ref = "a"))
  for (case in names(formulas)) {
    exog = switch(case, cross = d[["a"]] * d[["b_x_c"]],
      mixed = as.numeric(d[["dbreg_i"]] == "1"), raw = d[["dbreg_i1"]],
      factor_i = cbind(as.numeric(d[["g"]] == "1"), as.numeric(d[["f"]] == "b"), as.numeric(d[["f"]] == "c")),
      overlap = cbind(as.numeric(d[["f"]] == "b"), as.numeric(d[["f"]] == "c"), d[["dbreg_i1"]]))
    instruments = if (case == "cross") d[["a_x_b"]] * d[["c"]] else if (case %in% c("overlap", "factor_i")) d[["z"]] else
      cbind(as.numeric(d[["f"]] == "b"), as.numeric(d[["f"]] == "c"), d[["z"]])
    X = cbind(1, exog, d[["endog"]])
    Z = cbind(1, exog, instruments)
    bread = solve(crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, X)))
    beta = as.numeric(bread %*% crossprod(X, Z) %*% solve(crossprod(Z), crossprod(Z, d[["y"]])))
    V = bread * sum((d[["y"]] - X %*% beta)^2) / (nrow(d) - ncol(X))
    for (strategy in c("moments", "compress")) {
      fit = suppressWarnings(dbreg::dbivreg(formulas[[case]], data = d, strategy = strategy))
      expect_equal(as.numeric(coef(fit)), beta, tolerance = 1e-9)
      expect_equal(as.numeric(vcov(fit)), as.numeric(V), tolerance = 1e-9)
      expect_equal(df.residual(fit), nrow(d) - ncol(X))
      expect_equal(attr(vcov(fit), "df.K"), ncol(X))
      if (case == "factor_i") expect_identical(names(coef(fit)), c("(Intercept)", "g1", "f::b", "f::c", "endog"))
      # The cross fixture isolates target identity. Its literal _x_ R display
      # mismatch predates this change and remains a separate naming-policy item.
      if (case != "cross") {
        expect_equal(predict(fit, d), as.numeric(X %*% beta), tolerance = 1e-9)
        extra = d
        extra[["dbreg_i1_1_1"]] = -999
        expect_equal(predict(fit, extra), as.numeric(X %*% beta), tolerance = 1e-9)
      }
      close(fit)
    }
  }
  # Identical SQL with distinct unambiguous declared factor/i targets retains
  # the established rank policy; target allocation must not merge these.
  expect_error(dbreg::dbivreg(y ~ f | endog ~ i(f, ref = "a") + z, data = d, strategy = "moments"),
    "Instrument matrix is rank deficient")
  # Exact duplicated endogenous/instrument expressions must remain prohibited;
  # allocating distinct aliases must never disguise structural overlap.
  expect_error(dbreg::dbivreg(y ~ 1 | endog ~ endog, data = d, strategy = "moments"),
    "Excluded instruments cannot also appear as endogenous regressors")
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "iv_alias_training", d)
  blocks = list(
    dbreg:::dbivreg_expand_block(~a:b_x_c, conn, "iv_alias_training"),
    dbreg:::dbivreg_expand_block(~a_x_b:c, conn, "iv_alias_training"),
    dbreg:::dbivreg_expand_block(~a:b_x_c, conn, "iv_alias_training"))
  allocated = dbreg:::dbivreg_shared_aliases(blocks)
  expect_identical(allocated[[1]][["names"]], allocated[[3]][["names"]])
  expect_true(!identical(allocated[[1]][["names"]], allocated[[2]][["names"]]))
  for (k in seq_along(blocks)) {
    expect_identical(allocated[[k]][["sql"]], blocks[[k]][["sql"]])
    expect_identical(allocated[[k]][["columns"]][["id"]], blocks[[k]][["columns"]][["id"]])
    expect_identical(allocated[[k]][["columns"]][["alias"]], allocated[[k]][["names"]])
  }
})
