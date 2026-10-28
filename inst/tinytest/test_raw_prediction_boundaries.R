# Independent raw-prediction boundaries; raw labels never define the oracle X.
local({
  check = function(name, current, expected, tolerance = 1e-8) {
    expect_equal(current, expected, tolerance = tolerance, info = name)
  }
  run = function(mode) {
    set.seed(62481)
    n = if (mode == "composite") 144L else 120L
    d = data.frame(id = seq_len(n), x = rnorm(n), w = 1 + (seq_len(n) %% 5L) / 4)
    if (mode == "labels") {
      d[["grp"]] = rep(c("a", "b"), length.out = n)
      d[["grpb"]] = rnorm(n)
      d[["g"]] = rep(seq_len(4L), each = n / 4L)
      G = diag(4L)[d[["g"]], , drop = FALSE]
      X = cbind(1, as.numeric(d[["grp"]] == "b"), d[["grpb"]], G[, -1L, drop = FALSE])
      d[["y"]] = as.numeric(X %*% c(.2, 1.7, -.6, .3, -.4, .8)) + rnorm(n, sd = .15)
      formula = y ~ grp + grpb | g
      strategy = "compress"
      slopes = 2:3
    } else {
      keys = c(1, 1 + 2 * .Machine[["double.eps"]], 2)
      d[["g"]] = rep(keys, length.out = n)
      code = match(d[["g"]], keys)
      d[["x"]] = d[["x"]] + .4 * code
      if (mode == "composite") {
        d[["h"]] = rep(rep(1:2, each = 3L), length.out = n)
        code = code + 3L * (d[["h"]] - 1L)
        G = diag(6L)[code, , drop = FALSE]
        formula = y ~ x | g ^ h
        strategy = "compress"
      } else {
        G = diag(3L)[code, , drop = FALSE]
        formula = y ~ x | g
        strategy = mode
      }
      d[["y"]] = .6 * d[["x"]] + (code %% 3L) - .4 * code + rnorm(n, sd = .15)
      if (mode == "mundlak") {
        wx = rowsum(d[["w"]] * d[["x"]], code, reorder = TRUE)[, 1L]
        ww = rowsum(d[["w"]], code, reorder = TRUE)[, 1L]
        X = cbind(1, d[["x"]], (wx / ww)[code])
        slopes = seq_len(3L)
      } else {
        X = cbind(1, d[["x"]], G[, -1L, drop = FALSE])
        slopes = 2L
      }
    }
    reference = lm.wfit(X, d[["y"]], d[["w"]])
    V = solve(crossprod(X, d[["w"]] * X)) * sum(d[["w"]] * reference[["residuals"]]^2) / (n - ncol(X))
    conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    tables = sort(DBI::dbListTables(conn))
    model = dbreg(formula, data = d, conn = conn, strategy = strategy, weights = "w", retain = "none", vcov = "iid")
    on.exit(close(model), add = TRUE)
    check("beta", as.numeric(coef(model, fe = mode == "mundlak")), as.numeric(reference[["coefficients"]][slopes]))
    check("V", as.numeric(vcov(model, fe = mode == "mundlak")), as.numeric(V[slopes, slopes, drop = FALSE]))
    check("K", as.numeric(attr(vcov(model), "df.K")), as.numeric(ncol(X)))
    check("df", as.numeric(df.residual(model)), as.numeric(n - ncol(X)))
    if (mode %in% c("composite", "labels")) {
      cells = model[["data"]]
      if (mode == "labels") {
        cellG = diag(4L)[match(as.character(cells[["g"]]), as.character(1:4)), , drop = FALSE]
        cellX = cbind(1, as.numeric(cells[["grp"]] == "b"), cells[["grpb"]], cellG[, -1L, drop = FALSE])
      } else {
        cellcode = match(cells[["g"]], keys) + 3L * (cells[["h"]] - 1L)
        cellG = diag(6L)[cellcode, , drop = FALSE]
        cellX = cbind(1, cells[["x"]], cellG[, -1L, drop = FALSE])
      }
      check("stored_points", predict(model), as.numeric(cellX %*% reference[["coefficients"]]))
    }
    check("raw_R_points", predict(model, d), as.numeric(reference[["fitted.values"]]))
    close(model)
    check("caller_valid", DBI::dbIsValid(conn), TRUE)
    check("caller_tables", sort(DBI::dbListTables(conn)), tables)
    invisible(NULL)
  }
  for (mode in c("composite", "labels", "demean", "mundlak")) run(mode)
})
