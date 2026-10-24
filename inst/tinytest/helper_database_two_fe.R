dbreg_test_database_two_fe = function(conn) {
  before = dbglm_validation_tables(conn)
  original = getOption("dbreg.prediction_max_entries")
  on.exit(options(dbreg.prediction_max_entries = original), add = TRUE)
  set.seed(561)
  n = 960L
  d = data.frame(id = seq_len(n), g = rep(1:32, each = 30),
    h = sample(1:8, n, TRUE), x = rnorm(n), z = rnorm(n),
    w = rep(c(.25, .5, 1), length.out = n), o = .1*sin(seq_len(n)))
  d[["h"]][d[["g"]] > 16] = d[["h"]][d[["g"]] > 16] + 8
  d[["h"]][d[["g"]] == 1] = rep(1:2, 15)
  d[["tag"]] = paste0("group:'\\ café ", d[["g"]])
  d[["twox"]] = 2*d[["x"]] + sin(d[["g"]]) + cos(d[["h"]])
  d[["xs"]] = 1e6*d[["x"]]
  d[["y"]] = 2 + .3*d[["x"]] - .2*d[["z"]] + sin(d[["g"]]) + cos(d[["h"]]) + rnorm(n, sd = .1)
  d[["b"]] = as.numeric(d[["id"]] %% 3 == 0)
  d[["count"]] = as.numeric(1 + d[["id"]] %% 4)
  table = dbglm_validation_upload(conn, d)
  models = new.env(parent = emptyenv())
  on.exit({
    for (model in as.list(models)) close(model)
    dbglm_validation_drop(conn, table)
  }, add = TRUE)
  rows = unique(c(seq.int(1L, n, by = 30L), 2L, 13L, n))
  newdata = d[rows, setdiff(names(d), c("y", "b", "count", "w")), drop = FALSE]
  # g=1 never occurs with h=8, but that combination is estimable. h=16
  # belongs to the other disconnected component and must remain unidentified.
  novel = newdata[1:3, , drop = FALSE]
  novel[["g"]] = 1L
  novel[["tag"]] = d[["tag"]][[1L]]
  novel[["h"]] = c(8L, 16L, 100L)
  inventory = dbglm_validation_tables(conn)
  check = function(glm = FALSE, family = NULL, strategy = "demean") {
    response = if (!glm) "y" else if (family[["family"]] == "poisson") "count" else "b"
    message("[Two-FE SQL prediction] ", if (glm) family[["family"]] else strategy)
    fml = stats::as.formula(paste0(response, "~x+z|g^tag+h"), env = baseenv())
    args = list(fml = fml, conn = conn, table = table, weights = "w", retain = "prediction")
    if (glm) {
      args[["family"]] = family
      args[["offset"]] = "o"
      args[["fixef_rm"]] = "none"
    } else args[["strategy"]] = strategy
    options(dbreg.prediction_max_entries = original)
    portable = do.call(if (glm) dbglm else dbreg, args)
    options(dbreg.prediction_max_entries = 64)
    fit = do.call(if (glm) dbglm else dbreg, args)
    models[[as.character(length(models) + 1L)]] = fit
    expect_identical(fit[["fitted_state"]][["resources"]][["storage"]], "database")
    expect_equal(length(fit[["fitted_state"]][["resources"]][["tables"]]), 2L)
    expect_equal(length(setdiff(dbglm_validation_tables(conn), inventory)), 2L)
    expect_equal(coef(fit), coef(portable), tolerance = 1e-8)
    expect_equal(vcov(fit), vcov(portable), tolerance = 1e-8)
    expect_equal(nobs(fit), nobs(portable))
    prediction = predict(fit, newdata)
    expect_equal(prediction, predict(portable, newdata), tolerance = 1e-7)
    reference = if (glm) fixest::feglm(fml, d, weights = ~w, offset = ~o,
      family = family, fixef.rm = "none", notes = FALSE) else
        fixest::feols(fml, d, weights = ~w, fixef.rm = "none", notes = FALSE)
    expect_equal(prediction, as.numeric(predict(reference, newdata, type = "link")), tolerance = 1e-6)
    selected = c(7L, 2L, 11L)
    expect_equal(predict(fit, newdata[selected, , drop = FALSE]), prediction[selected], tolerance = 1e-10)
    novel_prediction = predict(fit, novel)
    expect_true(is.finite(novel_prediction[[1L]]) && all(is.na(novel_prediction[2:3])))
    expect_equal(novel_prediction[[1L]], predict(portable, novel)[[1L]], tolerance = 1e-7)
    restored = unserialize(serialize(fit, NULL))
    expect_error(predict(restored, newdata), "live owning session")
    close(restored)
    expect_equal(predict(fit, newdata), prediction, tolerance = 1e-10)
    close(fit)
    close(portable)
    expect_error(predict(fit, newdata), "released")
    expect_identical(dbglm_validation_tables(conn), inventory)
  }
  check(strategy = "demean")
  check(strategy = "compress")
  check(glm = TRUE, family = stats::binomial("logit"))
  check(glm = TRUE, family = stats::poisson())

  # FE-only GLMs need neither slope columns nor a within-X matrix. Exercise
  # both the diagonal and two-block SQL dictionaries with zero-column designs.
  for (specification in list(list(fml = b ~ 1 | g^tag, family = stats::binomial("logit")),
    list(fml = count ~ 1 | g^tag + h, family = stats::poisson()))) {
    options(dbreg.prediction_max_entries = 64)
    fit = dbglm(specification[["fml"]], conn = conn, table = table, weights = "w", offset = "o",
      family = specification[["family"]], retain = "prediction", fixef_rm = "none")
    models[[as.character(length(models) + 1L)]] = fit
    reference = fixest::feglm(specification[["fml"]], d, weights = ~w, offset = ~o,
      family = specification[["family"]], fixef.rm = "none", notes = FALSE)
    expect_identical(fit[["fitted_state"]][["resources"]][["storage"]], "database")
    expect_equal(predict(fit, newdata), as.numeric(predict(reference, newdata, type = "link")), tolerance = 1e-6)
    close(fit)
    expect_identical(dbglm_validation_tables(conn), inventory)
  }

  options(dbreg.prediction_max_entries = 64)
  fit = dbreg(y ~ x + twox + z | g^tag + h, conn = conn, table = table,
    weights = "w", strategy = "demean", retain = "prediction")
  models[[as.character(length(models) + 1L)]] = fit
  expect_true(all(is.finite(predict(fit, newdata))))
  changed = newdata[1:3, , drop = FALSE]
  changed[["twox"]][2L] = changed[["twox"]][2L] + 1
  prediction = predict(fit, changed)
  expect_true(is.finite(prediction[[1L]]) && is.na(prediction[[2L]]) && is.finite(prediction[[3L]]))
  # A new combination must still satisfy the absorbed regressor relation.
  novel_relation = novel[1L, , drop = FALSE]
  novel_relation[["twox"]] = 2*novel_relation[["x"]] + sin(novel_relation[["g"]]) + cos(novel_relation[["h"]])
  expect_true(is.finite(predict(fit, novel_relation)))
  novel_relation[["twox"]] = novel_relation[["twox"]] + 1
  expect_true(is.na(predict(fit, novel_relation)))
  expected = predict(fit, newdata)
  scaled = dbreg(y ~ xs + twox + z | g^tag + h, conn = conn, table = table,
    weights = "w", strategy = "demean", retain = "prediction")
  models[[as.character(length(models) + 1L)]] = scaled
  expect_equal(predict(scaled, newdata), expected, tolerance = 1e-7)
  expect_true(is.na(predict(scaled, changed)[[2L]]))
  close(scaled)
  dbglm_validation_drop(conn, table)
  newdata[["y"]] = -999
  expect_equal(predict(fit, newdata), expected, tolerance = 1e-10)
  close(fit)
  expect_identical(dbglm_validation_tables(conn), before)

  # A retention failure must not leave either dictionaries or row intermediates.
  table = dbglm_validation_upload(conn, d)
  inventory = dbglm_validation_tables(conn)
  expect_error(dbreg(y ~ x + z | g^tag + h, conn = conn, table = table, weights = "w",
    strategy = "compress", retain = "prediction", fe_max_iter = 1L), "recovery did not converge")
  expect_identical(dbglm_validation_tables(conn), inventory)
  expect_true(DBI::dbIsValid(conn))
}
