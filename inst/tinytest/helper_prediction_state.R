# Both engines run the same learned-effect and estimability checks. Prediction
# inputs omit the outcome and fitting weights, and include zero-weight rows.
dbglm_test_prediction_state = function(conn) {
  before = dbglm_validation_tables(conn)
  set.seed(709)
  n = 600L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n),
    g = rep(1:8, length.out = n), h = rep(1:5, each = 120),
    k = rep(1:3, each = 200), t = rnorm(n), cat = rep(c("a", "b", "c"), 200),
    w = rep(c(.25, .5, 1), 200), off = rnorm(n, sd = .08))
  d[["w"]][c(1, 17)] = 0
  eta = -.2 + .3 * d[["x"]] - .1 * d[["z"]] + .1 * sin(d[["g"]]) + .1 * cos(d[["h"]]) + d[["off"]]
  d[["yb"]] = rbinom(n, 1, plogis(eta))
  d[["yp"]] = rpois(n, exp(eta))
  source = dbglm_validation_upload(conn, d)
  tables = source
  on.exit(for (table in rev(tables)) dbglm_validation_drop(conn, table))
  inventory = sort(c(before, source))
  check = function(fml, family, method = "auto") {
    fit = dbglm(fml, conn = conn, table = source, family = family, weights = "w", offset = "off",
      fixef_rm = "none", retain = "prediction", control = dbglm_control(projection_method = method))
    ref = suppressWarnings(dbglm_reference_feglm(fml, d, family = family, weights = ~w, offset = ~off, fixef.rm = "none"))
    newdata = d[rev(seq_len(n)), setdiff(names(d), c("yb", "yp", "w")), drop = FALSE]
    found = predict(fit, newdata)
    expected = unname(predict(ref, newdata, type = "link"))
    expect_equal(found, expected, tolerance = 1e-6, info = paste(deparse(fml), family[["link"]], method))
    expect_equal(predict(fit, newdata, type = "response"), family[["linkinv"]](expected), tolerance = 1e-6)
    expect_equal(predict(fit, newdata[1, , drop = FALSE]), found[1], tolerance = 1e-10)
    expect_equal(predict(fit, newdata[c(3, 1, 3), , drop = FALSE]), found[c(3, 1, 3)], tolerance = 1e-10)
    expect_equal(predict(fit, newdata[FALSE, , drop = FALSE]), numeric())
    expect_equal(predict(unserialize(serialize(fit, NULL)), newdata), found, tolerance = 1e-12)
    # No environments, connections or functions may hide in mathematical state.
    values_only = function(value) {
      if (is.environment(value) || is.function(value) || typeof(value) == "externalptr") return(FALSE)
      if (is.list(value)) return(all(vapply(value, values_only, logical(1))))
      TRUE
    }
    expect_true(values_only(fit[["fitted_state"]]))
    expect_identical(dbglm_validation_tables(conn), inventory)
    expect_error(predict(fit), "newdata")
    fit
  }
  for (family in list(binomial("logit"), binomial("cloglog"), binomial("probit"), poisson())) {
    y = if (family[["family"]] == "poisson") "yp" else "yb"
    fit = check(as.formula(paste(y, "~ x + z + factor(cat) + I(x^2) | g + h")), family)
  }
  check(yb ~ x + z | g + h, binomial("cloglog"), method = "sql")
  check(yp ~ 1 | g, poisson())
  check(yp ~ x | g[t] + h, poisson())
  check(yb ~ x | g + h + k, binomial("logit"))
  check(yp ~ x | g^h, poisson())
  # Saved state is independent of later source mutation or a live connection.
  saved = fit
  input = d[1:4, setdiff(names(d), c("yb", "yp", "w")), drop = FALSE]
  expected = predict(saved, input)
  dbglm_validation_drop(conn, source)
  tables = setdiff(tables, source)
  expect_equal(predict(saved, input), expected, tolerance = 1e-12)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(dbglm_validation_tables(conn), before)

  # Known marginal levels do not identify every cross-component combination.
  cdata = data.frame(g = rep(1:4, each = 80), h = rep(c(1, 2, 3, 4), each = 80),
    x = rep(c(-1, 0, 1, 2), 80), y = rep(c(0, 1, 2, 1), 80))
  source = dbglm_validation_upload(conn, cdata)
  tables = c(tables, source)
  fit = dbglm(y ~ x | g + h, conn = conn, table = source, family = poisson(), retain = "prediction", fixef_rm = "none")
  known = cdata[c(1, 81), c("g", "h", "x")]
  expect_true(all(is.finite(predict(fit, known))))
  crossed = known
  crossed[["h"]] = rev(crossed[["h"]])
  expect_true(all(is.na(predict(fit, crossed))))
  crossed[["g"]][1] = 999
  expect_true(all(is.na(predict(fit, crossed))))
  known[["x"]][2] = NA_real_
  expect_true(is.finite(predict(fit, known)[1]) && is.na(predict(fit, known)[2]))
  expect_error(predict(fit, known["x"]), "missing columns")
  # A regressor relation holding only in training must also hold for a new row.
  cdata[["z"]] = 2 * cdata[["x"]]
  source2 = dbglm_validation_upload(conn, cdata)
  tables = c(tables, source2)
  collinear = dbglm(y ~ x + z | g, conn = conn, table = source2, family = poisson(), retain = "prediction", fixef_rm = "none")
  rows = cdata[1:2, c("g", "x", "z")]
  expect_true(all(is.finite(predict(collinear, rows))))
  rows[["z"]][2] = rows[["z"]][2] + 1
  expect_true(is.finite(predict(collinear, rows)[1]) && is.na(predict(collinear, rows)[2]))
  # A bounded-state failure must release private objects and preserve the source.
  original = getOption("dbreg.prediction_max_entries")
  on.exit(options(dbreg.prediction_max_entries = original), add = TRUE)
  inventory = dbglm_validation_tables(conn)
  options(dbreg.prediction_max_entries = 1)
  expect_error(dbglm(y ~ x | g, conn = conn, table = source, family = poisson(), retain = "prediction", fixef_rm = "none"), "prediction state exceeds")
  options(dbreg.prediction_max_entries = original)
  expect_identical(dbglm_validation_tables(conn), inventory)
}

dbglm_test_prediction_keys = function(conn) {
  before = dbglm_validation_tables(conn)
  d = data.frame(g = rep(1:4, each = 40), x = rep(c(-1, 0, 1, 2), 40),
    y = rep(c(0, 1, 2, 1), 40))
  d[["y"]] = d[["y"]] + d[["g"]] - 1
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("prediction_keys_"))
  on.exit({dbglm_validation_drop(conn, table); dbglm_validation_drop(conn, source)})
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  sql = paste0("SELECT *, CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS big,",
    "CAST(1.0 AS DOUBLE)+g*CAST(2.220446049250313e-16 AS DOUBLE) AS close_key,",
    "CASE WHEN g IN (1,2) THEN 'a|b' ELSE 'a' END AS part1,",
    "CASE WHEN g IN (1,3) THEN '' ELSE '|b' END AS part2 FROM ", q(source))
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table, sql, c(names(d), "big", "close_key", "part1", "part2")))
  newdata = d[c(1, 41, 81, 121), c("g", "x")]
  newdata[["big"]] = bit64::as.integer64("9007199254740992") + bit64::as.integer64(newdata[["g"]])
  newdata[["close_key"]] = 1 + newdata[["g"]] * 2.220446049250313e-16
  newdata[["part1"]] = ifelse(newdata[["g"]] %in% c(1, 2), "a|b", "a")
  newdata[["part2"]] = ifelse(newdata[["g"]] %in% c(1, 3), "", "|b")
  for (fml in list(y ~ x | big, y ~ x | close_key, y ~ x | part1^part2)) {
    fit = dbglm(fml, conn = conn, table = table, family = poisson(), retain = "prediction", fixef_rm = "none")
    expected = predict(fixest::feglm(y ~ x | g, d, family = poisson(), fixef.rm = "none", notes = FALSE), newdata, type = "link")
    expect_equal(predict(fit, newdata), unname(expected), tolerance = 1e-6)
    expect_identical(length(fit[["fitted_state"]][["specification"]][["dictionaries"]][[1]][["key"]]), 4L)
  }
  expect_identical(dbglm_validation_tables(conn), sort(c(before, source, table)))
}

dbglm_test_prediction_weight_scale = function(conn) {
  d = data.frame(g = rep(1:2, each = 20),
    h = rep(1:2, 20),
    y = c(rep(0, 16), rep(1, 4), rep(0, 4), rep(1, 16)),
    w = rep(c(1e-20, 1), each = 20))
  source = dbglm_validation_upload(conn, d)
  on.exit(dbglm_validation_drop(conn, source))
  for (family in list(binomial("logit"), binomial("cloglog"), binomial("probit"), poisson())) {
    for (fml in list(y ~ 1 | g, y ~ 1 | g + h)) {
      fit = dbglm(fml, conn = conn, table = source, family = family, weights = "w", retain = "prediction")
      expect_equal(predict(fit, data.frame(g = 1:2, h = 1), type = "response"), c(.2, .8), tolerance = 1e-7)
      expect_true(fit[["converged"]] && fit[["score_max"]] <= fit[["control"]][["score_tol"]])
    }
  }
}
