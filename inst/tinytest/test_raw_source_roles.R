# Raw source identities are opaque, including strings resembling interactions.
local({
  set.seed(29813)
  n = 96L
  d = data.frame(id = seq_len(n), x = rnorm(n), w = rep(c(.5, 1, 2), length.out = n))
  d[["A B"]] = d[["x"]]
  d[["y"]] = .3 + .7 * d[["x"]] + rnorm(n, sd = .4)
  raw_clusters = c("cluster value", "cluster^value", "cluster:value", "cluster\\value", "cluster`value", "A_B_x_C")
  for (j in seq_along(raw_clusters)) d[[raw_clusters[[j]]]] = as.integer((seq_len(n) - 1L) %/% (j + 1L) %% 8L + 1L)
  d[["future cluster"]] = sample.int(9L, n, TRUE)
  request = function(column) stats::as.formula(call("~", as.name(column)), env = baseenv())
  oracle = function(X, y, w, cluster) {
    fit = lm.wfit(X, y, w)
    bread = solve(crossprod(X, w * X))
    score = rowsum(X * (w * fit[["residuals"]]), cluster, reorder = FALSE)
    meat = crossprod(score)
    bare = bread %*% meat %*% bread
    G = nrow(score)
    correction = G / (G - 1) * (nrow(X) - 1) / (nrow(X) - ncol(X))
    list(beta = fit[["coefficients"]], fitted = fit[["fitted.values"]],
      residual = fit[["residuals"]], V = bare * correction, bare = bare)
  }
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "raw_role_rows", d)
  original_tables = sort(DBI::dbListTables(conn))
  collect_rows = function(output) {
    on.exit(close(output), add = TRUE)
    values = as.data.frame(output)
    column = output[["prediction_column"]]
    if (is.null(column)) column = output[["value_column"]]
    if (is.null(column)) column = "id"
    values[[column]][match(d[["id"]], values[["id"]])]
  }
  X = cbind(1, d[["x"]])
  for (column in raw_clusters) for (kind in c("OLS", "Gaussian")) {
    for (retention in if (column %in% c("cluster`value", "A_B_x_C")) c("inference", "all") else "inference") {
      vc = request(column)
      fit = if (kind == "OLS") dbreg(y ~ `A B`, conn = conn, table = "raw_role_rows",
        weights = "w", vcov = vc, retain = retention, row_id = "id", strategy = "moments") else
        dbglm(y ~ `A B`, conn = conn, table = "raw_role_rows", family = gaussian(),
          weights = "w", vcov = vc, retain = retention, row_id = "id", strategy = "moments")
      ref = oracle(X, d[["y"]], d[["w"]], d[[column]])
      expect_equal(unname(coef(fit)), unname(ref[["beta"]]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit)), as.numeric(ref[["V"]]), tolerance = 1e-8)
      expect_identical(fit[["inference_result"]][["clusters"]], list(column))
      expect_identical(fit[["cluster_var"]], column)
      stored = fit[["fitted_state"]][["inference"]][["default_vcov"]]
      expect_identical(all.vars(stored), column)
      expect_identical(environment(stored), baseenv())
      expect_equal(as.numeric(vcov(fit, ssc = list(K.adj = FALSE, G.adj = FALSE))),
        as.numeric(ref[["bare"]]), tolerance = 1e-8)
      future = oracle(X, d[["y"]], d[["w"]], d[["future cluster"]])
      expect_equal(as.numeric(vcov(fit, vcov = ~`future cluster`)), as.numeric(future[["V"]]), tolerance = 1e-8)
      if (kind == "Gaussian") {
        expect_identical(all.vars(fit[["call"]][["vcov"]]), column)
        expect_identical(environment(fit[["call"]][["vcov"]]), baseenv())
      }
      if (retention == "all") expect_equal(as.numeric(collect_rows(fitted(fit))), as.numeric(ref[["fitted"]]), tolerance = 1e-8)
      close(fit)
      expect_true(DBI::dbIsValid(conn))
      expect_identical(sort(DBI::dbListTables(conn)), original_tables)
    }
  }

  # Outcome/weight/offset names must address the immutable original rows.
  d[["A_B_x_Y"]] = d[["y"]]
  d[["A_B_x_W"]] = d[["w"]]
  d[["A_B_x_O"]] = rep(c(-.1, .2), length.out = n)
  d[["A_B_x_G"]] = rep(seq_len(4L), length.out = n)
  d[["A_B_x_H"]] = as.integer((seq_len(n) - 1L) %/% 4L %% 3L + 1L)
  d[["A_B_x_T"]] = rnorm(n)
  d[["A_B_x_X"]] = runif(n)
  DBI::dbWriteTable(conn, "raw_role_rows", d, overwrite = TRUE)
  for (strategy in c("moments", "compress")) {
    fit = dbreg(A_B_x_Y ~ `A B`, conn = conn, table = "raw_role_rows",
      weights = "A_B_x_W", offset = "A_B_x_O", vcov = ~A_B_x_C,
      retain = "all", row_id = "id", strategy = strategy, diagnostics = TRUE)
    ref = oracle(X, d[["A_B_x_Y"]] - d[["A_B_x_O"]], d[["A_B_x_W"]], d[["A_B_x_C"]])
    expect_equal(unname(coef(fit)), unname(ref[["beta"]]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(ref[["V"]]), tolerance = 1e-8)
    expect_equal(as.numeric(collect_rows(fitted(fit))), as.numeric(ref[["fitted"]] + d[["A_B_x_O"]]), tolerance = 1e-8)
    expect_equal(as.numeric(collect_rows(residuals(fit))), as.numeric(ref[["residual"]]), tolerance = 1e-8)
    expect_identical(fit[["yvar"]], "A_B_x_Y")
    expect_identical(fit[["weights"]], "A_B_x_W")
    expect_identical(fit[["cluster_var"]], "A_B_x_C")
    expect_equal(collect_rows(case.names(fit)), seq_len(n))
    detail = fit[["sample_diagnostics"]][["detail"]]
    expect_true(all(c("A_B_x_Y", "A B", "A_B_x_W", "A_B_x_C", "A_B_x_O") %in% detail[["variable"]]))
    expect_true(all(detail[["n_missing"]] == 0))
    close(fit)
  }

  # FE components and trend operations retain raw names and original row order.
  G = model.matrix(~factor(d[["A_B_x_G"]]) - 1)
  joint = cbind(1, d[["x"]], G[, -1, drop = FALSE], G * d[["A_B_x_T"]])
  y = d[["A_B_x_Y"]] - d[["A_B_x_O"]]
  ref = oracle(joint, y, d[["A_B_x_W"]], d[["A_B_x_C"]])
  for (strategy in c("demean", "compress")) {
    fit = dbreg(A_B_x_Y ~ `A B` | A_B_x_G[A_B_x_T], conn = conn, table = "raw_role_rows",
      weights = "A_B_x_W", offset = "A_B_x_O", vcov = ~A_B_x_C, retain = "all", row_id = "id", strategy = strategy)
    expect_equal(unname(coef(fit)), unname(ref[["beta"]][2L]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(ref[["V"]][2L, 2L]), tolerance = 1e-7)
    expect_equal(predict(fit, d), as.numeric(ref[["fitted"]] + d[["A_B_x_O"]]), tolerance = 1e-7)
    expect_equal(as.numeric(collect_rows(fitted(fit))), as.numeric(ref[["fitted"]] + d[["A_B_x_O"]]), tolerance = 1e-7)
    expect_identical(fit[["varying_slopes"]][[1L]], "A_B_x_T")
    expect_identical(fit[["fitted_state"]][["specification"]][["dictionaries"]][[1L]][["columns"]], "A_B_x_G")
    expect_equal(collect_rows(case.names(fit)), seq_len(n))
    close(fit)
  }

  D = model.matrix(~interaction(d[["A_B_x_G"]], d[["A_B_x_H"]], drop = TRUE))
  composite_X = cbind(D[, 1L, drop = FALSE], d[["x"]], D[, -1L, drop = FALSE])
  composite_ref = oracle(composite_X, y, d[["A_B_x_W"]], d[["A_B_x_C"]])
  composite = dbreg(A_B_x_Y ~ `A B` | A_B_x_G^A_B_x_H, conn = conn, table = "raw_role_rows",
    weights = "A_B_x_W", offset = "A_B_x_O", vcov = ~A_B_x_C,
    retain = "all", row_id = "id", strategy = "compress")
  expect_equal(as.numeric(coef(composite)), as.numeric(composite_ref[["beta"]][2L]), tolerance = 1e-7)
  expect_equal(as.numeric(vcov(composite)), as.numeric(composite_ref[["V"]][2L, 2L]), tolerance = 1e-7)
  expect_equal(predict(composite, d), as.numeric(composite_ref[["fitted"]] + d[["A_B_x_O"]]), tolerance = 1e-7)
  expect_identical(composite[["composite_fes"]][[1L]], c("A_B_x_G", "A_B_x_H"))
  expect_equal(as.character(composite[["fe"]]), names(composite[["composite_fes"]]))
  close(composite)
  cells = dbreg(A_B_x_Y ~ `A B` | A_B_x_G, conn = conn, table = "raw_role_rows",
    strategy = "compress", data_only = TRUE)
  expect_true(all(c("A B", "A_B_x_G") %in% names(cells)))
  expect_equal(sum(cells[["n"]]), n)
  expect_equal(sort(cells[["A B"]]), sort(d[["x"]]))
  expect_equal(sort(as.character(cells[["A_B_x_G"]])), sort(as.character(d[["A_B_x_G"]])))

  # An i() dictionary does not turn the same literal raw numeric source into
  # categorical main-effect coding merely because its name contains _x_.
  d[["A_B_x_P"]] = rep(c(10L, 20L, 40L), length.out = n)
  DBI::dbWriteTable(conn, "raw_role_rows", d, overwrite = TRUE)
  mixed = dbreg(A_B_x_Y ~ `A B` + A_B_x_P + i(A_B_x_P, ref = 10),
    conn = conn, table = "raw_role_rows", weights = "A_B_x_W", strategy = "compress")
  mixed_X = cbind(1, d[["x"]], d[["A_B_x_P"]],
    as.numeric(d[["A_B_x_P"]] == 20), as.numeric(d[["A_B_x_P"]] == 40))
  mixed_ref = lm.wfit(mixed_X, d[["A_B_x_Y"]], d[["w"]])
  expect_equal(predict(mixed, d), as.numeric(mixed_ref[["fitted.values"]]), tolerance = 1e-7)
  expect_true("A_B_x_P" %in% names(coef(mixed, collin = TRUE)))
  close(mixed)

  # IV offset outcome restoration consumes the same raw source identity.
  d[["z"]] = rnorm(n)
  d[["endog"]] = .7 * d[["z"]] + .2 * d[["x"]] + rnorm(n, sd = .2)
  d[["A_B_x_Y"]] = .4 + .5 * d[["x"]] + .8 * d[["endog"]] + d[["A_B_x_O"]] + rnorm(n, sd = .2)
  DBI::dbWriteTable(conn, "raw_role_rows", d, overwrite = TRUE)
  C = cbind(1, d[["x"]], d[["endog"]])
  Z = cbind(1, d[["x"]], d[["z"]])
  WZ = d[["w"]] * Z
  projection = crossprod(C, WZ) %*% solve(crossprod(Z, WZ))
  expected = solve(projection %*% crossprod(Z, d[["w"]] * C),
    projection %*% crossprod(Z, d[["w"]] * (d[["A_B_x_Y"]] - d[["A_B_x_O"]])))
  iv = dbivreg(A_B_x_Y ~ `A B` | endog ~ z, conn = conn, table = "raw_role_rows",
    weights = "A_B_x_W", offset = "A_B_x_O", strategy = "moments")
  expect_equal(as.numeric(coef(iv)), as.numeric(expected), tolerance = 1e-8)
  expect_identical(iv[["yvar"]], "A_B_x_Y")
  close(iv)

  # Bins source identities/inner completed clusters use whole aliases, not tokens.
  breaks = seq(0, 1, length.out = 5L)
  bin = cut(d[["A_B_x_X"]], breaks, include.lowest = TRUE, right = FALSE)
  B = model.matrix(~bin)
  for (with_fe in c(FALSE, TRUE)) {
    fml = if (with_fe) A_B_x_Y ~ A_B_x_X + `A B` | A_B_x_G else A_B_x_Y ~ A_B_x_X + `A B`
    bins = dbbinsreg(fml, data = d, nbins = 4L, binspos = breaks, points = c(0, 0),
      strategy = "compress", sample_fit = FALSE, vcov = ~A_B_x_C, plot = FALSE)
    XX = cbind(B, d[["x"]], if (with_fe) G[, -1, drop = FALSE])
    binref = oracle(XX, d[["A_B_x_Y"]], rep(1, n), d[["A_B_x_C"]])
    model = bins[["model"]]
    expect_equal(as.numeric(coef(model, fe = with_fe)), as.numeric(binref[["beta"]]), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(model, fe = with_fe)), as.numeric(binref[["V"]]), tolerance = 1e-7)
    expect_identical(bins[["opt"]][["x_var"]], "A_B_x_X")
    expect_identical(bins[["opt"]][["y_var"]], "A_B_x_Y")
    expect_identical(model[["inference_result"]][["clusters"]], list("A_B_x_C"))
    if (with_fe) expect_identical(model[["fe_labels"]], "A_B_x_G") else {
      at = cbind(1, diag(4L)[, -1, drop = FALSE], 0)
      expect_equal(bins[["points"]][["fit"]], as.numeric(at %*% binref[["beta"]]), tolerance = 1e-7)
    }
    close(model)
  }
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), original_tables)
})
