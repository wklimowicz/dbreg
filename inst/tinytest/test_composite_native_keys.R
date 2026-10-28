# Typed tuple controls do not infer native types from value text or labels.
local({
  adjacent = data.frame(g = c(1, 1 + 2 * .Machine[["double.eps"]]), h = c("a|b", "a|b"))
  expect_equal(length(unique(dbreg:::tuple_key(adjacent, c("g", "h"), c("numeric", "character")))), 2L)
  text = data.frame(a = c("a|b", "a"), b = c("c", "b|c"))
  expect_equal(length(unique(dbreg:::tuple_key(text, c("a", "b"), c("character", "character")))), 2L)
  values = data.frame(day = as.Date(c("2020-01-01", "2020-01-02", NA)), bit = c(FALSE, TRUE, TRUE))
  keys = dbreg:::tuple_key(values, c("day", "bit"), c("date", "logical"))
  expect_equal(length(unique(keys[1:2])), 2L)
  expect_true(is.na(keys[[3L]]))
  integer = data.frame(g = c("9007199254740993", "9007199254740994"), h = c(1L, 1L))
  expect_equal(length(unique(dbreg:::tuple_key(integer, c("g", "h"), c("integer", "integer")))), 2L)
  expect_error(dbreg:::tuple_key(data.frame(g = c(2^53 + 2, 2^53 + 4)), "g", "integer"), "must be exact")
  expect_error(dbreg:::tuple_key(adjacent, c("g", "h"), c(NA_character_, "character")), "complete")
})

# A BIGINT composite component also used as numeric RHS must keep its original
# native key before the structural DOUBLE projection. Exercise numeric DBI mode.
if (requireNamespace("bit64", quietly = TRUE)) local({
  set.seed(62482)
  n = 120L
  labels = c("9007199254740992", "9007199254740993", "9007199254740994")
  d = data.frame(id = seq_len(n), g = rep(labels, length.out = n),
    h = rep(rep(1:2, each = 3L), length.out = n), x = rnorm(n), w = 1 + (seq_len(n) %% 5L) / 4)
  code = match(d[["g"]], labels) + 3L * (d[["h"]] - 1L)
  G = diag(6L)[code, , drop = FALSE]
  X = cbind(1, d[["x"]], G[, -1L, drop = FALSE])
  d[["y"]] = .4 * d[["x"]] + code / 3 + rnorm(n, sd = .1)
  reference = lm.wfit(X, d[["y"]], d[["w"]])
  V = solve(crossprod(X, d[["w"]] * X)) * sum(d[["w"]] * reference[["residuals"]]^2) / (n - ncol(X))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")), bigint = "numeric")
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(conn, "native_components", d)
  DBI::dbExecute(conn, "CREATE TEMP VIEW native_component_input AS SELECT id,CAST(g AS BIGINT) AS g,h,x,w,y FROM native_components")
  tables = sort(DBI::dbListTables(conn))
  model = dbreg(y ~ x + g | g ^ h, conn = conn, table = "native_component_input", weights = "w",
    strategy = "compress", retain = "none", vcov = "iid")
  on.exit(close(model), add = TRUE)
  expect_identical(names(coef(model)), "x")
  expect_equal(as.numeric(coef(model)), as.numeric(reference[["coefficients"]][2L]), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(model)), as.numeric(V[2L, 2L]), tolerance = 1e-8)
  expect_equal(attr(vcov(model), "df.K"), ncol(X))
  expect_true(is.na(coef(model, collin = TRUE)[["g"]]))
  lookup = model[["composite_levels"]][[1L]]
  expect_equal(length(lookup), 6L)
  expect_identical(attr(lookup, "types"), c("integer", "integer"))
  expect_false(anyDuplicated(names(model[["data"]])) > 0L)
  expect_true(is.character(model[["data"]][["g"]]))
  expect_equal(length(unique(model[["data"]][["g"]])), 3L)
  expect_false(any(grepl("^dbreg_composite_type", names(model[["data"]]))))
  newdata = d[c("x", "g", "h")]
  newdata[["g"]] = bit64::as.integer64(newdata[["g"]])
  expect_equal(predict(model, newdata), as.numeric(reference[["fitted.values"]]), tolerance = 1e-8)
  cells = model[["data"]]
  cellcode = match(cells[["g"]], labels) + 3L * (cells[["h"]] - 1L)
  cellX = cbind(1, cells[["x"]], diag(6L)[cellcode, -1L, drop = FALSE])
  expect_equal(predict(model), as.numeric(cellX %*% reference[["coefficients"]]), tolerance = 1e-8)
  expect_identical(dbreg:::tuple_key(newdata, c("g", "h"), attr(lookup, "types")),
    dbreg:::tuple_key(d, c("g", "h"), attr(lookup, "types")))
  lossy = newdata
  lossy[["g"]] = as.numeric(lossy[["g"]])
  expect_error(predict(model, lossy), "must be exact")
  close(model)
  cells_only = dbreg(y ~ x + g | g ^ h, conn = conn, table = "native_component_input", weights = "w",
    strategy = "compress", data_only = TRUE)
  expect_false(any(grepl("^dbreg_compress_rhs", names(cells_only))))
  expect_equal(sum(cells_only[["n"]]), n)
  expect_equal(length(unique(cells_only[["g"]])), 3L)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), tables)
})
