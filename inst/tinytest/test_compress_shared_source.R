# One raw source can supply distinct ordinary and FE computational roles.
local({
  set.seed(84217)
  n = 180L
  d = data.frame(id = seq_len(n), g = rep(c(10L, 20L, 40L), length.out = n),
    x = rnorm(n), z = rnorm(n), w = rep(c(.5, 1, 2), length.out = n),
    cl = rep(seq_len(12L), each = 15L))
  d[["dbreg_compress_fe1__"]] = seq_len(n) + 300L
  d[["g value"]] = d[["g"]]
  d[["y"]] = .7 * d[["x"]] + .1 * d[["g"]] + rnorm(n, sd = .2)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "shared_roles", d)
  before = sort(DBI::dbListTables(conn))
  collect = function(output) {
    on.exit(close(output), add = TRUE)
    rows = as.data.frame(output)
    column = output[["prediction_column"]]
    if (is.null(column)) column = output[["value_column"]]
    as.numeric(rows[[column]][match(d[["id"]], rows[["id"]])])
  }
  check = function(model, X, slopes) {
    ref = lm.wfit(X, d[["y"]], d[["w"]])
    bread = solve(crossprod(X, d[["w"]] * X))
    scores = rowsum(X * (d[["w"]] * ref[["residuals"]]), d[["cl"]], reorder = FALSE)
    V = bread %*% crossprod(scores) %*% bread * 12 / 11 * (n - 1) / (n - ncol(X))
    expect_equal(as.numeric(coef(model)), as.numeric(ref[["coefficients"]][slopes]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(model)), as.numeric(V[slopes, slopes, drop = FALSE]), tolerance = 1e-8)
    expect_equal(attr(vcov(model), "df.K"), ncol(X))
    expect_equal(predict(model, d), as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    expect_equal(collect(predict(model, conn = conn, table = "shared_roles", row_id = "id")),
      as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    expect_equal(collect(fitted(model)), as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    invisible(ref)
  }
  for (raw in c("g", "g value")) {
    rhs = call("+", as.name("x"), as.name(raw))
    fml = stats::as.formula(call("~", as.name("y"), call("|", rhs, as.name(raw))))
    model = dbreg(fml, conn = conn, table = "shared_roles", weights = "w",
      strategy = "compress", vcov = ~cl, retain = "all", row_id = "id")
    X = model.matrix(~x + factor(g), d)
    ref = check(model, X, 2L)
    layout = model[["linear_design"]]
    private = names(layout[["fe_sources"]])
    expect_equal(length(private), 1L)
    expect_identical(unname(layout[["fe_sources"]]), raw)
    expect_identical(unname(layout[["source_map"]][private]), raw)
    expect_identical(layout[["factor_roles"]][[private]], "fe")
    expect_true(private %in% names(model[["factor_levels"]]))
    expect_equal(length(model[["factor_levels"]][[private]]), 3L)
    expect_true(is.factor(model[["data"]][[private]]))
    expect_true(is.numeric(model[["data"]][[raw]]))
    expect_false(anyDuplicated(names(model[["data"]])) > 0L)
    ordinary = dbreg:::model_vars(layout[["rhs"]])
    ordinary = ordinary[unname(layout[["source_map"]][ordinary]) == raw]
    positions = which(vapply(layout[["columns"]][["recipe"]], function(recipe)
      identical(recipe[["kind"]], "numeric") && identical(recipe[["variable"]], ordinary), logical(1)))
    expect_equal(length(positions), 1L)
    expect_identical(layout[["columns"]][["role"]][positions], "slope")
    expect_false(layout[["columns"]][["id"]][positions] %in% layout[["fitted_ids"]])
    raw_label = if (raw == "g") "g" else "`g value`"
    expect_true(is.na(coef(model, collin = TRUE)[[raw_label]]))
    minimal = d[c("x", raw)]
    expect_equal(predict(model, minimal), as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    minimal[[private]] = rev(d[["g"]])
    expect_equal(predict(model, minimal), as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    expect_error(predict(model, minimal[setdiff(names(minimal), raw)]), "missing")
    expect_equal(predict(model, model[["data"]]),
      as.numeric(coef(model, fe = TRUE)[["(Intercept)"]]) +
        as.numeric(model[["data"]][["x"]]) * as.numeric(coef(model)[["x"]]) +
        as.numeric(model.matrix(~factor(model[["data"]][[private]]))[, -1L, drop = FALSE] %*%
          tail(as.numeric(coef(model, fe = TRUE)), 2L)), tolerance = 1e-8)
    close(model)
  }
  cells = dbreg(y ~ x + g | g, conn = conn, table = "shared_roles", strategy = "compress", data_only = TRUE)
  expect_false(any(grepl("^dbreg_compress_fe", names(cells))))
  expect_true(all(c("g", "x", "n") %in% names(cells)))
  expect_equal(sum(cells[["n"]]), n)

  # Categorical logical FE context determines the interaction's reference basis.
  d[["g"]] = rep(c("a", "b", "c"), length.out = n)
  DBI::dbWriteTable(conn, "shared_roles", d, overwrite = TRUE)
  categorical = dbreg(y ~ x + g + g:z | g, conn = conn, table = "shared_roles",
    weights = "w", vcov = ~cl, strategy = "compress", retain = "all", row_id = "id")
  Z = model.matrix(~g - 1, d)
  X = cbind(model.matrix(~x + g, d), Z[, -1L, drop = FALSE] * d[["z"]])
  ref = check(categorical, X, c(2L, 5L, 6L))
  expect_equal(length(coef(categorical)), 3L)
  close(categorical)

  # Private centered trend fields are rebuilt from original raw numeric input.
  d[["g"]] = rep(c(10L, 20L, 40L), length.out = n)
  d[["trend"]] = 2000 + rnorm(n)
  DBI::dbWriteTable(conn, "shared_roles", d, overwrite = TRUE)
  trend = dbreg(y ~ x + g | g[trend], conn = conn, table = "shared_roles",
    weights = "w", vcov = ~cl, strategy = "compress", retain = "all", row_id = "id")
  Z = model.matrix(~factor(g), d)
  X = cbind(Z[, 1L, drop = FALSE], d[["x"]], Z[, -1L, drop = FALSE],
    model.matrix(~factor(g) - 1, d) * (d[["trend"]] - 2000))
  check(trend, X, 2L)
  close(trend)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), before)
})
