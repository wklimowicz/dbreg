local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  old = options(dbreg.ap_solver = "auto", dbreg.ap_method = NULL)
  on.exit(options(old), add = TRUE)
  set.seed(11122)
  n = 600L
  d = data.frame(id = seq_len(n), x = sample(1:5, n, TRUE), z = rnorm(n),
    g = sample(c("a^b", "a", "", "é|", "NULL"), n, TRUE),
    h = sample(c("c", "b^c", "", "|/"), n, TRUE), k = sample(1:4, n, TRUE),
    w = runif(n, .2, 2), o = rnorm(n))
  d[["endog"]] = .8 * d[["z"]] + rnorm(n)
  # An independent integer Cartesian index avoids any delimiter oracle.
  d[["joint"]] = match(d[["g"]], unique(d[["g"]])) + length(unique(d[["g"]])) * (match(d[["h"]], unique(d[["h"]])) - 1)
  d[["y"]] = d[["o"]] + .7 * d[["x"]] + .5 * d[["endog"]] + sin(d[["joint"]]) + cos(d[["k"]]) + rnorm(n)
  DBI::dbWriteTable(conn, "joint_fixture", d)
  initial = sort(DBI::dbListTables(conn))
  adjusted = d
  adjusted[["y"]] = d[["y"]] - d[["o"]]
  for (two in c(FALSE, TRUE)) for (strategy in c("demean", "compress")) for (kind in list("iid", "hc1", ~g)) {
    fml = if (two) y ~ x | g^h + k else y ~ x | g^h
    fit = dbreg(fml, conn = conn, table = "joint_fixture", weights = ~w, offset = ~o,
      strategy = strategy, vcov = kind, diagnostics = TRUE)
    ref = fixest::feols(if (two) y ~ x | joint + k else y ~ x | joint, adjusted, weights = ~w,
      vcov = if (identical(kind, "hc1")) "hetero" else kind,
      ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
    expect_identical(fit[["fe"]], if (two) c("g^h", "k") else "g^h")
    expect_equal(nobs(fit), n)
    expect_equal(fit[["sample_diagnostics"]][["n_fitted"]], n)
    expect_identical(fit[["sample_diagnostics"]][["detail"]][["variable"]][fit[["sample_diagnostics"]][["detail"]][["role"]] == "fe"], if (two) c("g", "h", "k") else c("g", "h"))
    # Full prediction uses the training tuple coding, even when newdata is
    # permuted or only contains a subset of its levels.
    pred = predict(fit, d)
    D = model.matrix(if (two) ~x + factor(joint) + factor(k) else ~x + factor(joint), d)
    reference = lm.wfit(D, adjusted[["y"]], d[["w"]])
    expect_equal(as.numeric(pred), as.numeric(D %*% reference[["coefficients"]] + d[["o"]]), tolerance = 1e-7)
    if (strategy == "compress") {
      take = seq.int(n, 1L, by = -7L)
      expect_equal(as.numeric(predict(fit, d[take, , drop = FALSE])), as.numeric(pred[take]), tolerance = 1e-8)
      restored = unserialize(serialize(fit, NULL))
      expect_equal(as.numeric(predict(restored, d[take, , drop = FALSE])), as.numeric(pred[take]), tolerance = 1e-8)
      unknown = d[1:3, , drop = FALSE]
      unknown[["g"]][2] = "unseen tuple"
      expect_identical(is.na(predict(fit, unknown)), c(FALSE, TRUE, FALSE))
    }
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  for (two in c(FALSE, TRUE)) for (kind in list("iid", "hc1", ~g)) {
    fit = dbivreg(if (two) y ~ x | g^h + k | endog ~ z else y ~ x | g^h | endog ~ z,
      conn = conn, table = "joint_fixture", weights = ~w, offset = ~o, vcov = kind)
    ref = fixest::feols(if (two) y ~ x | joint + k | endog ~ z else y ~ x | joint | endog ~ z,
      adjusted, weights = ~w, vcov = if (identical(kind, "hc1")) "hetero" else kind,
      ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
    ix = match(names(coef(fit)), sub("^fit_", "", names(coef(ref))))
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)[ix]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[ix, ix, drop = FALSE]), tolerance = 1e-8)
    expect_identical(fit[["fe"]], if (two) c("g^h", "k") else "g^h")
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # Native BIGINT tuples retain adjacent values above 2^53 without a numeric
  # key round trip. Original row_id exports retain those exact typed keys.
  DBI::dbExecute(conn, "CREATE TABLE big_joint AS SELECT * EXCLUDE(k),CAST('9007199254740992' AS BIGINT)+CAST(k AS BIGINT) AS k FROM joint_fixture")
  for (strategy in c("demean", "compress")) {
    fit = dbreg(y ~ x | g^k + h, conn = conn, table = "big_joint", weights = ~w,
      offset = ~o, strategy = strategy, row_id = c("id", "k"), sample_table = "joint_keys")
    ref = fixest::feols(y ~ x | g^k + h, adjusted, weights = ~w, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    keys = DBI::dbGetQuery(conn, "SELECT CAST(k AS VARCHAR) AS k FROM joint_keys ORDER BY id")[["k"]]
    expect_identical(keys, paste0("900719925474099", d[["k"]] + 2L))
    DBI::dbRemoveTable(conn, "joint_keys")
    exact = d
    exact[["k"]] = keys
    canonical = predict(fit, exact)
    exact[["k"]] = paste0("0", keys)
    expect_equal(as.numeric(predict(fit, exact)), as.numeric(canonical), tolerance = 1e-8)
    if (strategy == "compress") {
      raw = d
      raw[["k"]] = keys
      expect_true(all(is.finite(predict(fit, raw))))
      raw[["k"]] = as.numeric(keys)
      expect_error(predict(fit, raw), "exact|integer64|decimal strings")
    }
  }
  DBI::dbRemoveTable(conn, "big_joint")
  # Component nulls must be filtered, never converted into a valid joint ID.
  missing = d
  missing[["g"]][1:4] = NA_character_
  missing[["h"]][5:7] = NA_character_
  fit = dbreg(y ~ x | g^h + k, data = missing, weights = ~w, diagnostics = TRUE)
  expect_equal(nobs(fit), n - 7)
  expect_equal(fit[["sample_diagnostics"]][["detail"]][["n_missing"]][fit[["sample_diagnostics"]][["detail"]][["variable"]] %in% c("g", "h")], c(4, 3))
  # Cohort Mundlak means group actual tuples independently of training ranks.
  fit = dbreg(y ~ x | g^h + k, data = d, weights = ~w, offset = ~o, strategy = "mundlak")
  mean_by = function(group) ave(seq_len(n), group, FUN = function(i) weighted.mean(d[["x"]][i], d[["w"]][i]))
  X = cbind(1, d[["x"]], mean_by(d[["joint"]]), mean_by(d[["k"]]))
  ref = lm.wfit(X, adjusted[["y"]], d[["w"]])
  expect_equal(as.numeric(predict(fit, d)), as.numeric(X %*% ref[["coefficients"]] + d[["o"]]), tolerance = 1e-8)
  # A literal source name that resembles a tuple is a separate producer.
  d[["g^h"]] = as.numeric(d[["x"]])
  fit = dbreg(y ~ `g^h` | g^h + k, data = d, weights = ~w, strategy = "compress")
  expect_identical(names(coef(fit)), "`g^h`")
  expect_true(all(is.finite(predict(fit, d))))
  # Dictionary output aliases and paste() argument names cannot capture raw
  # source roles named id or collapse.
  roles = d
  roles[["id"]] = roles[["g"]]
  roles[["collapse"]] = roles[["h"]]
  fit = dbreg(y ~ x | id^collapse + k, data = roles, weights = ~w, strategy = "compress")
  expect_true(all(is.finite(predict(fit, roles))))
  for (limit in list(list(dbreg.projection_max_rows = 5), list(dbreg.projection_max_cells = 10), list(dbreg.projection_max_bytes = 1024))) {
    prior = options(limit)
    expect_error(dbreg(y ~ x | g^h + k, conn = conn, table = "joint_fixture", weights = ~w, strategy = "compress"), "prediction dictionaries exceed")
    options(prior)
    expect_identical(sort(DBI::dbListTables(conn)), initial)
  }
  # The two bounded dictionaries share one cumulative budget.
  prior = options(dbreg.projection_max_rows = 25)
  expect_error(dbreg(y ~ x | g^h + g^k, conn = conn, table = "joint_fixture", weights = ~w, strategy = "compress"), "prediction dictionaries exceed")
  options(prior)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  sql = capture.output({plan = dbreg(y ~ x | g^h + k, conn = conn, table = "joint_fixture", weights = ~w, sql_only = TRUE, strategy = "demean")})
  expect_true(inherits(plan, "dbreg_query_plan"))
  expect_identical(sort(DBI::dbListTables(conn)), initial)
  adjacent = data.frame(g = c(1, 1 + 2 * .Machine[["double.eps"]]), h = "a|b")
  expect_equal(length(unique(dbreg:::composite_key(adjacent, c("g", "h"), c("numeric", "character")))), 2L)
  dated = data.frame(day = as.Date(c("2020-01-01", "2020-01-02", NA)), bit = c(FALSE, TRUE, TRUE))
  keys = dbreg:::composite_key(dated, c("day", "bit"), c("date", "logical"))
  expect_equal(length(unique(keys[1:2])), 2L)
  expect_true(is.na(keys[3L]))
  expect_error(dbreg:::composite_key(adjacent, c("g", "h"), c(NA_character_, "character")), "complete")
  # Native keys shared with a numeric structural term keep distinct roles.
  labels = c("9007199254740992", "9007199254740993", "9007199254740994")
  m = 120L
  native = data.frame(g = rep(labels, length.out = m), h = rep(rep(1:2, each = 3L), length.out = m), x = rnorm(m), w = 1 + (seq_len(m) %% 5L) / 4)
  code = match(native[["g"]], labels) + 3L * (native[["h"]] - 1L)
  X = cbind(1, native[["x"]], diag(6L)[code, -1L, drop = FALSE])
  native[["y"]] = .4 * native[["x"]] + code / 3 + rnorm(m, sd = .1)
  ref = lm.wfit(X, native[["y"]], native[["w"]])
  V = solve(crossprod(X, native[["w"]] * X)) * sum(native[["w"]] * ref[["residuals"]]^2) / (m - ncol(X))
  DBI::dbWriteTable(conn, "native_components", native)
  DBI::dbExecute(conn, "CREATE TEMP VIEW native_component_input AS SELECT CAST(g AS BIGINT) AS g,h,x,w,y FROM native_components")
  for (kind in list("iid", "hc1", ~h)) {
    fit = dbreg(y ~ x + g | g^h, conn = conn, table = "native_component_input", weights = ~w, strategy = "compress", vcov = kind)
    expect_identical(names(coef(fit)), "x")
    expect_equal(as.numeric(coef(fit)), as.numeric(ref[["coefficients"]][2L]), tolerance = 1e-8)
    if (identical(kind, "iid")) expect_equal(as.numeric(vcov(fit)), V[2L, 2L], tolerance = 1e-8)
    expect_true("g" %in% fit[["collin.var"]])
    expect_true(is.character(fit[["data"]][["g"]]))
    expect_equal(length(unique(fit[["data"]][["g"]])), 3L)
    expect_false(anyDuplicated(names(fit[["data"]])) > 0L)
    expect_equal(as.numeric(predict(fit, native)), as.numeric(ref[["fitted.values"]]), tolerance = 1e-8)
    expect_true(all(is.finite(predict(fit))))
  }
  cells = dbreg(y ~ x + g | g^h, conn = conn, table = "native_component_input", weights = ~w, strategy = "compress", data_only = TRUE)
  expect_equal(sum(cells[["n"]]), m)
  expect_equal(length(unique(cells[["g"]])), 3L)
  expect_false(any(grepl("dbreg_joint_rhs", names(cells))))
  DBI::dbExecute(conn, "DROP VIEW native_component_input")
  DBI::dbRemoveTable(conn, "native_components")
  # No oversized raw key is fetched before preflight, including data_only.
  wide = data.frame(y = rnorm(20), x = rnorm(20), g = rep(strrep("a", 65535), 20), h = rep(1:2, 10))
  DBI::dbWriteTable(conn, "wide_components", wide)
  state = new.env(parent = emptyenv());state[["fetches"]] = 0L
  execute = dbreg:::execute_compress_strategy
  environment(execute) = new.env(parent = environment(execute))
  environment(execute)[["read_compressed_moments"]] = function(...) {
    state[["fetches"]] = state[["fetches"]] + 1L
    dbreg:::read_compressed_moments(...)
  }
  api = dbreg
  environment(api) = new.env(parent = environment(api))
  environment(api)[["execute_compress_strategy"]] = execute
  for (data_only in c(FALSE, TRUE)) {
    expect_error(api(y ~ x | g^h, conn = conn, table = "wide_components", strategy = "compress", data_only = data_only), "key strings exceed")
    expect_identical(state[["fetches"]], 0L)
    expect_identical(sort(DBI::dbListTables(conn)), sort(c(initial, "wide_components")))
  }
  DBI::dbRemoveTable(conn, "wide_components")
  # Canonical count/moment names cannot replace a raw tuple component.
  roles = d
  roles[["n"]] = roles[["k"]]
  roles[["sum_w"]] = roles[["h"]]
  for (kind in list("iid", ~g)) {
    fit = dbreg(y ~ x | n^sum_w + g, data = roles, weights = ~w, strategy = "compress", vcov = kind)
    ref = fixest::feols(y ~ x | k^h + g, roles, weights = ~w, vcov = kind, ssc = fixest::ssc(K.fixef = "full"), fixef.rm = "none", notes = FALSE)
    expect_equal(nobs(fit), n)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
    expect_false(anyDuplicated(names(fit[["data"]])) > 0L)
    expect_equal(sum(fit[["data"]][["n"]]), n)
    expect_true(all(c("n", "sum_w") %in% names(attr(fit[["data"]], "dbreg_source_columns"))))
    expect_true(all(is.finite(predict(fit, roles))))
    expect_equal(as.numeric(predict(fit, fit[["data"]])), as.numeric(predict(fit)), tolerance = 1e-8)
  }
  cells = dbreg(y ~ x | n^sum_w + g, data = roles, weights = ~w, strategy = "compress", data_only = TRUE)
  expect_equal(sum(cells[["n"]]), n)
  expect_false(anyDuplicated(names(cells)) > 0L)
  expect_equal(length(unique(cells[[attr(cells, "dbreg_source_columns")[["n"]]]])), 4L)
  for (fml in list(y ~ x | g:h, y ~ x | g*h, y ~ x | log(g), y ~ x | g^2, y ~ x | g[h])) {
    expect_error(dbreg(fml, data = d), "fixed effects")
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
