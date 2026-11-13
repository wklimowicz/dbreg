# Independent reference-coded spline/control/FE design proves plotted levels,
# their entire covariance contrasts, and the returned model's full joint V.
local({
  main = function() {
    conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
    n = 160L
    index = seq_len(n)
    data = data.frame(x = (index - .5) / n,
      z = 2 + sin(index * .71) + cos(index * .37),
      g = rep(sprintf("g%02d", seq_len(8L)), length.out = n),
      cluster = rep(seq_len(20L), length.out = n))
    breaks = seq(0, 1, length.out = 5L)
    bin = cut(data[["x"]], breaks, include.lowest = TRUE, right = FALSE)
    knots = as.numeric(tapply(data[["x"]], bin, max))[-4L]
    x_mean = as.numeric(tapply(data[["x"]], bin, mean))
    spline = cbind(x1_spline = data[["x"]], vapply(knots,
      function(knot) pmax(data[["x"]] - knot, 0), numeric(n)))
    colnames(spline) = c("x1_spline", sprintf("knot%d_pow1", seq_along(knots)))
    data[["y"]] = 4 + drop(spline %*% c(.7, -.4, .25, .15)) +
      .8 * data[["z"]] + .35 * (as.integer(factor(data[["g"]])) - 1L) +
      .12 * sin(index * 1.17)
    reference = function(fe, kind) {
      nuisance = if (fe) outer(data[["g"]], sprintf("g%02d", 2:8), `==`) * 1 else NULL
      if (fe) colnames(nuisance) = paste0("g", sprintf("g%02d", 2:8))
      X = cbind("(Intercept)" = 1, spline, z = data[["z"]], nuisance)
      fit = stats::lm.fit(X, data[["y"]])
      expect_identical(fit[["rank"]], ncol(X))
      K = ncol(X)
      bread = solve(crossprod(X))
      scores = X * fit[["residuals"]]
      G = length(unique(data[["cluster"]]))
      V = if (kind == "iid") bread * sum(fit[["residuals"]]^2) / (n - K) else if (kind == "hc1") {
        n / (n - K) * bread %*% crossprod(scores) %*% bread
      } else {
        grouped = rowsum(scores, data[["cluster"]], reorder = FALSE)
        G / (G - 1) * (n - 1) / (n - K) * bread %*% crossprod(grouped) %*% bread
      }
      Z = matrix(0, length(x_mean), ncol(X), dimnames = list(NULL, colnames(X)))
      Z[, "(Intercept)"] = 1
      Z[, "x1_spline"] = x_mean
      for (j in seq_along(knots)) Z[, sprintf("knot%d_pow1", j)] = pmax(x_mean - knots[[j]], 0)
      # z=0 and reference FE g01=0, matching the unconstrained level target.
      values = drop(Z %*% fit[["coefficients"]])
      standard_error = sqrt(diag(Z %*% V %*% t(Z)))
      df_t = if (kind == "cluster") G - 1 else n - K
      list(beta = fit[["coefficients"]], V = V, values = values, se = standard_error,
        lwr = values - stats::qt(.975, df_t) * standard_error,
        upr = values + stats::qt(.975, df_t) * standard_error, K = K, G = G,
        df = n - K, df_t = df_t)
    }
    check = function(fe, kind, ci = TRUE, nonstandard = FALSE, strategy = "auto") {
      input = data
      vars = c(y = "y", x = "x", z = "z", g = "g", cluster = "cluster")
      if (nonstandard) {
        renamed = c(y = "response value", x = "running value", z = "control value",
          g = "panel value", cluster = "cluster value")
        names(input)[match(vars, names(input))] = unname(renamed)
        vars = renamed
      }
      rhs = call("+", as.name(vars[["x"]]), as.name(vars[["z"]]))
      if (fe) rhs = call("|", rhs, as.name(vars[["g"]]))
      fml = stats::as.formula(call("~", as.name(vars[["y"]]), rhs), env = baseenv())
      request = if (kind == "cluster") stats::as.formula(call("~", as.name(vars[["cluster"]])), env = baseenv()) else kind
      DBI::dbWriteTable(conn, "bins_spline_input", input, overwrite = TRUE)
      before = sort(DBI::dbListTables(conn))
      captured = new.env(parent = emptyenv())
      captured[["warnings"]] = character()
      fit = withCallingHandlers(dbbinsreg(fml, conn = conn, table = "bins_spline_input",
        points = c(1L, 1L), nbins = 4L, binspos = breaks,
        ci = ci, cb = FALSE, vcov = request, sample_fit = FALSE,
        strategy = strategy, plot = FALSE), warning = function(warning) {
          message = conditionMessage(warning)
          captured[["warnings"]] = c(captured[["warnings"]], message)
        })
      expect_true(all(grepl("compression ineffective", captured[["warnings"]], fixed = TRUE)))
      expected = reference(fe, kind)
      if (nonstandard) {
        names(expected[["beta"]])[names(expected[["beta"]]) == "z"] = "`control value`"
        dimnames(expected[["V"]]) = rep(list(names(expected[["beta"]])), 2L)
      }
      model = fit[["model"]]
      beta = stats::coef(model, fe = TRUE)
      V = stats::vcov(model, fe = TRUE)
      expect_identical(names(beta), names(expected[["beta"]]))
      expect_identical(rownames(V), names(expected[["beta"]]))
      expect_identical(colnames(V), names(expected[["beta"]]))
      expect_equal(as.numeric(beta), as.numeric(expected[["beta"]]), tolerance = 1e-8)
      expect_equal(as.numeric(V), as.numeric(expected[["V"]]), tolerance = 1e-8)
      expect_equal(as.numeric(attr(V, "n_params")), expected[["K"]])
      expect_equal(as.numeric(attr(V, "df.t")), expected[["df_t"]])
      expect_equal(model[["df_residual"]], expected[["df"]])
      if (kind == "cluster") expect_equal(as.numeric(attr(V, "n_clusters")), expected[["G"]])
      expect_equal(fit[["points"]][["x"]], x_mean, tolerance = 1e-12)
      expect_equal(fit[["points"]][["fit"]], expected[["values"]], tolerance = 1e-8)
      if (ci) {
        expect_equal(fit[["points"]][["se"]], expected[["se"]], tolerance = 1e-8)
        expect_equal(fit[["points"]][["lwr"]], expected[["lwr"]], tolerance = 1e-8)
        expect_equal(fit[["points"]][["upr"]], expected[["upr"]], tolerance = 1e-8)
      } else expect_false(any(c("se", "lwr", "upr") %in% names(fit[["points"]])))
      if (fe || strategy == "compress") expect_identical(model[["strategy"]], "compress")
      expect_identical(fit[["opt"]][["ci"]], ci)
      expect_true(DBI::dbIsValid(conn))
      expect_identical(sort(DBI::dbListTables(conn)), before)
      invisible(fit)
    }
    for (kind in c("iid", "hc1", "cluster")) {
      check(FALSE, kind)
      check(TRUE, kind)
    }
    check(TRUE, "iid", ci = FALSE)
    # One control proves both the nonstandard source boundary and forwarding an
    # explicit no-FE compression request after row-level spline construction.
    renamed = check(FALSE, "iid", nonstandard = TRUE, strategy = "compress")
    expect_identical(renamed[["opt"]][["x_var"]], "running value")
    expect_identical(renamed[["opt"]][["y_var"]], "response value")
  }
  main()
})
