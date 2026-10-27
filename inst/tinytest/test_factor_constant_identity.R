# A real category named _constant cannot acquire the synthetic constant's role.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  for (level in c("_constant", "single")) for (forced in c(FALSE, TRUE)) for (intercept in c(FALSE, TRUE)) {
    d = data.frame(y = seq_len(12) + rep(c(-0.1, 0.2), 6), cat = rep(level, 12))
    DBI::dbWriteTable(conn, "constant_identity", d, overwrite = TRUE)
    formula = if (forced) {
      if (intercept) y ~ factor(cat) else y ~ 0 + factor(cat)
    } else if (intercept) y ~ cat else y ~ 0 + cat
    design = dbreg::sql_model_matrix(formula, conn, "constant_identity", intercept = intercept)
    expect_identical(length(design[["select_exprs"]]), 1L)
    expect_identical(length(design[["aliases"]]), 1L)
    expect_identical(length(design[["col_names"]]), 1L)
    expect_identical(nrow(design[["columns"]]), 1L)
    expect_identical(length(design[["columns"]][["recipe"]]), 1L)
    expect_identical(anyDuplicated(design[["columns"]][["id"]]), 0L)
    # The legacy public mapping retains both rows,without confusing their roles.
    variable = if (forced) "factor(cat)" else "cat"
    expect_identical(sum(design[["level_aliases"]][["var"]] == variable), 2L)
    query = paste0("SELECT ", design[["select_exprs"]], " AS ", design[["aliases"]], " FROM constant_identity")
    actual = DBI::dbGetQuery(conn, query)
    expect_identical(ncol(actual), 1L)
    expect_equal(as.numeric(actual[[1L]]), rep(1, 12))
    for (strategy in if (forced) c("moments", "compress") else "moments") {
      model = dbreg::dbreg(formula, data = d, conn = conn, strategy = strategy)
      expect_equal(as.numeric(coef(model)), mean(d[["y"]]), tolerance = 1e-10)
      expected_V = sum((d[["y"]] - mean(d[["y"]]))^2) / 11 / 12
      expect_equal(as.numeric(vcov(model)), expected_V, tolerance = 1e-10)
      expect_equal(attr(vcov(model), "df.K"), 1)
      expect_equal(df.residual(model), 11)
      expect_equal(predict(model, d[1:3, ]), rep(mean(d[["y"]]), 3), tolerance = 1e-10)
      close(model)
    }
  }
})
