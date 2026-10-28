# FE nuisance span takes precedence over structurally aliased regressors.
local({
  set.seed(28416)
  n = 240L
  d = expand.grid(g = seq_len(4L), h = seq_len(3L), period = seq_len(5L), replication = seq_len(4L))
  d[["id"]] = seq_len(n)
  d[["x"]] = rnorm(n)
  d[["w"]] = rep(c(.5, 1, 2), length.out = n)
  d[["cl"]] = rep(seq_len(12L), each = 20L)
  d[["x2"]] = d[["x"]] + .3 * d[["g"]]
  d[["y"]] = .7 * d[["x"]] + sin(d[["g"]]) + cos(d[["h"]]) + .2 * d[["period"]] + rnorm(n, sd = .3)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "absorbed_rank", d)
  original_tables = sort(DBI::dbListTables(conn))
  oracle = function(X, weights, kind, K, cluster = d[["cl"]]) {
    fit = lm.wfit(X, d[["y"]], weights)
    bread = solve(crossprod(X, weights * X))
    e = fit[["residuals"]]
    rss = sum(weights * e^2)
    V = if (kind == "iid") bread * rss / (n - K) else {
      score = X * (weights * e)
      if (kind == "cluster") score = rowsum(score, cluster, reorder = FALSE)
      meat = crossprod(score)
      correction = if (kind == "hc1") n / (n - K) else {
        G = length(unique(cluster))
        G / (G - 1) * (n - 1) / (n - K)
      }
      bread %*% meat %*% bread * correction
    }
    list(beta = fit[["coefficients"]], V = V, fitted = fit[["fitted.values"]], residual = e, rss = rss)
  }
  cases = list(
    list(formula = y ~ x + i(g) | g, reference = ~x + factor(g), slopes = 1L, weighted = TRUE),
    list(formula = y ~ x + g | g, reference = ~x + factor(g), slopes = 1L, weighted = FALSE),
    list(formula = y ~ x + factor(g) | g, reference = ~x + factor(g), slopes = 1L, weighted = TRUE),
    list(formula = y ~ x + i(g) | g + h, reference = ~x + factor(g) + factor(h), slopes = 1L, weighted = TRUE),
    list(formula = y ~ x + i(period) | g, reference = ~x + factor(period) + factor(g), slopes = 5L, weighted = FALSE)
  )
  collect = function(output) {
    on.exit(close(output), add = TRUE)
    rows = as.data.frame(output)
    column = output[["prediction_column"]]
    if (is.null(column)) column = output[["value_column"]]
    as.numeric(rows[[column]][match(d[["id"]], rows[["id"]])])
  }
  for (case_index in seq_along(cases)) {
    spec = cases[[case_index]]
    X = model.matrix(spec[["reference"]], d)
    # Reference matrices are explicitly full rank; K also equals declared FE policy
    # here because the primary two-way design is connected and fully crossed.
    expect_equal(qr(X)[["rank"]], ncol(X))
    slopes = seq_len(spec[["slopes"]]) + 1L
    w = if (spec[["weighted"]]) d[["w"]] else rep(1, n)
    for (kind in c("iid", "hc1", "cluster")) {
      model = dbreg(spec[["formula"]], conn = conn, table = "absorbed_rank",
        strategy = "compress", weights = if (spec[["weighted"]]) "w" else NULL,
        vcov = if (kind == "cluster") ~cl else kind,
        ssc = list(K.fixef = "full"), retain = if (case_index == 1L && kind == "iid") "all" else "none",
        row_id = "id")
      ref = oracle(X, w, kind, ncol(X))
      expect_equal(length(coef(model)), length(slopes))
      expect_equal(dim(vcov(model)), rep(length(slopes), 2L))
      expect_equal(as.numeric(coef(model)), as.numeric(ref[["beta"]][slopes]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(model)), as.numeric(ref[["V"]][slopes, slopes, drop = FALSE]), tolerance = 1e-8)
      expect_equal(attr(model[["vcov"]], "rss"), ref[["rss"]], tolerance = 1e-8)
      expect_equal(df.residual(model), n - ncol(X))
      expect_equal(attr(vcov(model), "df.K"), ncol(X))
      layout = model[["linear_design"]]
      expect_equal(length(layout[["slope_ids"]]), length(slopes))
      positions = match(layout[["fitted_ids"]], layout[["columns"]][["id"]])
      expect_identical(positions, sort(positions))
      expect_equal(predict(model, d), as.numeric(ref[["fitted"]]), tolerance = 1e-8)
      if (case_index < 5L) expect_true(anyNA(coef(model, collin = TRUE))) else
        expect_equal(names(coef(model)), c("x", paste0("period::", 2:5)))
      if (case_index == 1L && kind == "iid") {
        expect_equal(collect(fitted(model)), as.numeric(ref[["fitted"]]), tolerance = 1e-8)
        expect_equal(collect(residuals(model)), as.numeric(ref[["residual"]]), tolerance = 1e-8)
        sql_prediction = predict(model, conn = conn, table = "absorbed_rank", row_id = "id")
        expect_equal(collect(sql_prediction), as.numeric(ref[["fitted"]]), tolerance = 1e-8)
        # Two additional retained SSC actions, independently computed.
        none = vcov(model, ssc = list(K.fixef = "none"))
        none_ref = oracle(X, w, "iid", 1L)
        expect_equal(as.numeric(none), as.numeric(none_ref[["V"]][2L, 2L]), tolerance = 1e-8)
        expect_equal(attr(none, "df.K"), 1L)
        nested = vcov(model, vcov = ~g, ssc = list(K.fixef = "nonnested"))
        # Nonnested excludes the three FE contrasts, retaining the intercept;
        # unlike K.fixef="none", its declared K is one slope plus one intercept.
        nested_ref = oracle(X, w, "cluster", 2L, d[["g"]])
        expect_equal(as.numeric(nested), as.numeric(nested_ref[["V"]][2L, 2L]), tolerance = 1e-8)
        expect_equal(attr(nested, "df.K"), 2L)
        reduced = fixest::feols(y ~ x | g, data = d, weights = ~w, vcov = ~g,
          ssc = fixest::ssc(K.fixef = "nonnested"), fixef.rm = "none",
          fixef.tol = 1e-11, nthreads = 1)
        expect_equal(fixest::degrees_freedom(reduced, "k"), 2L)
        expect_equal(as.numeric(vcov(reduced)), as.numeric(nested_ref[["V"]][2L, 2L]), tolerance = 1e-8)
      }
      close(model)
    }
  }
  # Neither ordinary column is individually constant within g, but their
  # residual directions coincide. Preserve the first surviving compiler axis.
  X = model.matrix(~x + factor(g), d)
  ref = oracle(X, d[["w"]], "iid", ncol(X))
  partial = dbreg(y ~ x + x2 | g, conn = conn, table = "absorbed_rank",
    weights = "w", strategy = "compress", vcov = "iid")
  expect_identical(names(coef(partial)), "x")
  expect_equal(as.numeric(coef(partial)), as.numeric(ref[["beta"]][2L]), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(partial)), as.numeric(ref[["V"]][2L, 2L]), tolerance = 1e-8)
  expect_equal(predict(partial, d), as.numeric(ref[["fitted"]]), tolerance = 1e-8)
  expect_true(is.na(coef(partial, collin = TRUE)[["x2"]]))
  close(partial)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), original_tables)
})
