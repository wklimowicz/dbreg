# Cast numeric operands before products, including design interactions and
# outcome/weight powers; casting SUM(product) cannot prevent integer overflow.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")), bigint = "numeric")
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(conn, paste(
    "CREATE TABLE products AS SELECT CAST((id%11-5)*20000 AS INTEGER) AS z,",
    "CAST((id%11-5)*16000+(id%7-3)*12000 AS INTEGER) AS x,",
    "CAST((id%11-5)*10000+(id%13-6)*20000 AS INTEGER) AS y,",
    "CAST(id%3+1 AS INTEGER) AS w, id%9 AS cl,",
    "CASE WHEN id%2=0 THEN 'a' ELSE 'b' END AS f FROM range(180) t(id)"))
  d = DBI::dbGetQuery(conn, "SELECT * FROM products")
  # The reference must also promote R integers before formula multiplication.
  d[c("x", "z", "y", "w")] = lapply(d[c("x", "z", "y", "w")], as.numeric)
  for (type in c("INTEGER", "BIGINT", "DECIMAL(18,0)")) {
    table = paste0("(SELECT CAST(x AS ", type, ") AS x, CAST(z AS ", type,
      ") AS z, CAST(y AS ", type, ") AS y, w, cl, f FROM products) typed")
    # Test the SQL design's values directly against R, independently of fitting.
    design = sql_model_matrix(~x:z + f:x, conn, table)
    queried = tryCatch(DBI::dbGetQuery(conn, paste("SELECT", paste(design[["select_exprs"]], collapse = ","), "FROM", table)), error = function(e) e)
    expect_false(inherits(queried, "error"), info = type)
    expected = stats::model.matrix(~0+x:z+f:x, d)
    # f:x alone keeps both levels; match by the compiler's public column names.
    expected = expected[, c("x:z", "x:fa", "x:fb"), drop = FALSE]
    if (!inherits(queried, "error")) expect_equal(unname(as.matrix(queried)), unname(expected), tolerance = 1e-12)
    for (strategy in c("moments", "compress")) {
      for (weighted in c(FALSE, TRUE)) {
        for (kind in list("iid", "hc1", ~cl)) {
          for (iv in c(FALSE, TRUE)) {
            formula = if (iv) y ~ 0 | x ~ z else y ~ 0 + x:z
            estimator = if (iv) dbivreg else dbreg
            fit = tryCatch(estimator(formula, conn = conn, table = table, strategy = strategy,
              weights = if (weighted) "w" else NULL, vcov = kind), error = function(e) e)
            info = paste(type, strategy, weighted, deparse(kind), iv)
            expect_false(inherits(fit, "error"), info = paste(info, if (inherits(fit, "error")) conditionMessage(fit) else ""))
            if (inherits(fit, "error")) next
            ref = fixest::feols(formula, d, weights = if (weighted) ~w else NULL, vcov = kind)
            expect_equal(nobs(ref), 180)
            b = coef(ref); names(b) = sub("^fit_", "", names(b))
            # Compare in ordinary units so tiny interaction coefficients and
            # variances cannot pass merely by being near zero.
            units = if (iv) 1 else 1e5
            expect_equal(as.numeric(coef(fit))*units, as.numeric(b)*units, tolerance = 1e-7, info = info)
            expect_equal(as.numeric(vcov(fit))*units^2, as.numeric(vcov(ref))*units^2, tolerance = 1e-7, info = info)
          }
        }
      }
    }
  }
  # i() numeric-by operands must also be promoted before their moments.
  fit = tryCatch(dbreg(y ~ 0 + i(f, x), conn = conn, table = "products", strategy = "moments", vcov = "hc1"), error = function(e) e)
  expect_false(inherits(fit, "error"))
  if (!inherits(fit, "error")) {
    ref = fixest::feols(y ~ 0 + i(f, x), d, vcov = "hc1")
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-7)
  }
  expect_true(DBI::dbIsValid(conn))
  expect_identical(DBI::dbListTables(conn), "products")
})
