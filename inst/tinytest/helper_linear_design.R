# The same synthetic formula cases run on DuckDB and native Databricks.
dbreg_test_linear_backend = function(conn) {
  set.seed(862)
  n = 240L
  dat = data.frame(x = rnorm(n), a = rep(1:3, 80), b = rep(1:4, each = 60),
    g = sample(1:5, n, TRUE), grp = rep(c("control", "treated"), 120),
    cl = sample(1:12, n, TRUE), w = runif(n, 0.5, 2))
  dat[["y"]] = dat[["x"]] - 0.3 * dat[["x"]]^2 + dat[["a"]] + rnorm(n)
  dat[["id"]] = seq_len(n)
  source = dbglm_validation_upload(conn, dat)
  table = paste0(source, "_sample")
  on.exit({
    dbglm_validation_drop(conn, table)
    dbglm_validation_drop(conn, source)
  })
  # The coefficient-lookup uploader deliberately rejects missing numeric values.
  # Construct missing observations in SQL, as in the other native sample tests.
  columns = c("CASE WHEN id=1 THEN NULL ELSE x END AS x",
    "CASE WHEN id=2 THEN NULL ELSE b END AS b", setdiff(names(dat), c("x", "b")))
  dbreg:::create_temp_table_as(conn, table,
    paste("SELECT", paste(columns, collapse = ","), "FROM", source),
    dbreg:::detect_backend(conn)[["name"]])
  dat[["x"]][1] = NA_real_
  dat[["b"]][2] = NA_integer_
  reference_data = dat
  reference_data[["grp"]] = factor(dat[["grp"]], levels = c("treated", "control"))
  for (strategy in c("moments", "compress", "demean")) {
    fml = if (strategy == "moments") y ~ x + I(x^2) + grp else y ~ x + I(x^2) + grp | a ^ b + g
    fit = dbreg(fml, conn = conn, table = table, weights = "w", vcov = ~cl,
      strategy = strategy, factor_refs = c(grp = "treated"), diagnostics = TRUE, fe_tol = 1e-8)
    ref = fixest::feols(fml, reference_data, weights = ~w, vcov = ~cl,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE), notes = FALSE)
    names_ref = sub("^I\\(I\\((.*)\\)\\)$", "I(\\1)", names(coef(ref)))
    index = match(names(coef(fit)), names_ref)
    expect_equal(sort(names(coef(fit))), sort(names_ref))
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)[index]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[index, index]), tolerance = 1e-7)
    expect_equal(nobs(fit), nobs(ref))
    expect_equal(fit[["nobs_dropped"]], n - nobs(ref))
    # Match the estimation population for conditional demean predictions.
    sample = if (strategy == "moments") !is.na(dat[["x"]]) else complete.cases(dat)
    expect_equal(predict(fit, dat[sample, ]), as.numeric(predict(ref, reference_data[sample, ])), tolerance = 1e-6)
  }
}
