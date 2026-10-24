dbreg_test_database_prediction = function(conn) {
  before = dbglm_validation_tables(conn)
  original = getOption("dbreg.prediction_max_entries")
  on.exit(options(dbreg.prediction_max_entries = original), add = TRUE)
  n = 768L
  d = data.frame(id = seq_len(n), g = rep(seq_len(64), each = 12),
    x = sin(seq_len(n) * .71), z = cos(seq_len(n) * .31),
    w = rep(c(.25, .5, 1), length.out = n))
  d[["s"]] = ifelse(d[["g"]] %% 2 == 0, 1, 1 + 2^-52)
  d[["y"]] = 2 + .3*d[["x"]] - .2*d[["z"]] + sin(d[["g"]]) + .1*cos(d[["id"]])
  d[["b"]] = as.numeric(d[["id"]] %% 4 %in% c(1, 2))
  d[["count"]] = as.numeric(1 + d[["id"]] %% 4)
  d[["o"]] = .1 * sin(d[["id"]] * .4)
  d[["twox"]] = 2*d[["x"]]
  d[["tag"]] = paste0("grp:'\\ café ", d[["g"]])
  d[["day"]] = as.character(as.Date("2026-10-01") + d[["g"]] %% 5)
  d[["flag"]] = ifelse(d[["g"]] %% 2 == 0, "true", "false")
  source = dbglm_validation_upload(conn, d)
  table = basename(tempfile("database_prediction_"))
  retained = basename(tempfile("prediction_sample_"))
  models = new.env(parent = emptyenv())
  on.exit({
    for (model in as.list(models)) close(model)
    dbglm_validation_drop(conn, retained)
    dbglm_validation_drop(conn, table)
    dbglm_validation_drop(conn, source)
  }, add = TRUE)
  fields = names(d)
  fields[fields == "g"] = "CAST(9007199254740992 AS BIGINT)+CAST(g AS BIGINT) AS g"
  fields[fields == "w"] = "CASE WHEN id=2 THEN 0.0 ELSE w END AS w"
  fields[fields == "y"] = "CASE WHEN id=1 THEN CAST(NULL AS DOUBLE) ELSE y END AS y"
  fields[fields == "day"] = "CAST(day AS DATE) AS day"
  fields[fields == "flag"] = "CAST(flag AS BOOLEAN) AS flag"
  DBI::dbExecute(conn, dbreg:::sql_sample_table(conn, table,
    paste("SELECT", paste(fields, collapse = ","), "FROM", source), names(d)))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  newdata = DBI::dbGetQuery(conn, paste0("SELECT id,CAST(g AS STRING) AS g,s,x,z,o,twox,tag,day,flag FROM ", q(table), " ORDER BY id DESC"))
  # Cover all groups and excluded rows without repeating warehouse round trips
  # over every synthetic row under the deliberately tiny transfer budget.
  selected_rows = sort(unique(c(round(seq.int(1L, n, length.out = 96L)), n - 1L)))
  newdata = newdata[selected_rows, , drop = FALSE]
  reference_rows = d[as.numeric(newdata[["id"]]), , drop = FALSE]
  d[["y"]][1L] = NA_real_
  d[["w"]][2L] = 0
  inventory = dbglm_validation_tables(conn)
  check = function(glm = FALSE, family = NULL, strategy = "demean", export = FALSE) {
    response = if (!glm) "y" else if (family[["family"]] == "poisson") "count" else "b"
    fml = stats::as.formula(paste0(response, "~x+z|g^s"), env = baseenv())
    args = list(fml = fml, conn = conn, table = table, weights = "w", retain = "prediction")
    if (glm) {
      args[["family"]] = family
      args[["offset"]] = "o"
      args[["fixef_rm"]] = "none"
    } else args[["strategy"]] = strategy
    options(dbreg.prediction_max_entries = original)
    portable = do.call(if (glm) dbglm else dbreg, args)
    options(dbreg.prediction_max_entries = 64)
    if (export) {
      args[["row_id"]] = "id"
      args[["sample_table"]] = retained
    }
    fit = do.call(if (glm) dbglm else dbreg, args)
    models[[as.character(length(models) + 1L)]] = fit
    expect_identical(fit[["fitted_state"]][["resources"]][["storage"]], "database")
    expect_equal(coef(fit), coef(portable), tolerance = 1e-8)
    expect_equal(vcov(fit), vcov(portable), tolerance = 1e-8)
    expect_equal(nobs(fit), nobs(portable))
    prediction = predict(fit, newdata)
    expect_equal(prediction, predict(portable, newdata), tolerance = 1e-7)
    ref = if (glm) fixest::feglm(fml, d, weights = ~w, offset = ~o, family = family,
      fixef.rm = "none", notes = FALSE) else fixest::feols(fml, d, weights = ~w, fixef.rm = "none", notes = FALSE)
    expect_equal(prediction, as.numeric(predict(ref, reference_rows, type = "link")), tolerance = 1e-6)
    selected = newdata[c(7L, 1L, 24L), , drop = FALSE]
    expect_equal(predict(fit, selected), prediction[c(7L, 1L, 24L)], tolerance = 1e-10)
    selected[["g"]][2L] = "9007199254741999"
    selected[["x"]][3L] = NA_real_
    expect_true(is.finite(predict(fit, selected)[1L]) && all(is.na(predict(fit, selected)[2:3])))
    restored = unserialize(serialize(fit, NULL))
    expect_error(predict(restored, newdata), "live owning session")
    close(restored)
    expect_equal(predict(fit, newdata), prediction, tolerance = 1e-10)
    close(fit)
    expect_error(predict(fit, newdata), "released")
    close(fit)
    expect_true(DBI::dbIsValid(conn))
    expect_identical(dbglm_validation_tables(conn), sort(c(inventory, if (export) retained)))
    if (export) {
      expect_equal(as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", q(retained)))[["n"]]), nobs(fit))
      dbglm_validation_drop(conn, retained)
    }
    close(portable)
    fit
  }
  check(strategy = "demean", export = TRUE)
  check(strategy = "compress")
  check(glm = TRUE, family = stats::binomial("logit"))
  check(glm = TRUE, family = stats::poisson())

  # Typed joins also preserve arbitrary text, dates and logical components.
  fit = dbreg(y ~ x + z | g^s^tag^day^flag, conn = conn, table = table, weights = "w",
    strategy = "demean", retain = "prediction")
  models[[as.character(length(models) + 1L)]] = fit
  reference = fixest::feols(y ~ x + z | g, d, weights = ~w, fixef.rm = "none", notes = FALSE)
  expect_equal(predict(fit, newdata), as.numeric(predict(reference, reference_rows)), tolerance = 1e-6)
  unknown = newdata[1:3, , drop = FALSE]
  unknown[["tag"]][1L] = "new:'\\ key"
  unknown[["day"]] = as.character(unknown[["day"]])
  unknown[["day"]][2L] = "invalid-date"
  unknown[["flag"]][3L] = !unknown[["flag"]][3L]
  expect_true(all(is.na(predict(fit, unknown))))
  close(fit)

  # A regressor relation absorbed by the FE must still hold in newdata.
  options(dbreg.prediction_max_entries = 64)
  fit = dbreg(y ~ x + twox | g, conn = conn, table = table, weights = "w",
    strategy = "demean", retain = "prediction")
  models[[as.character(length(models) + 1L)]] = fit
  expect_true(all(is.finite(predict(fit, newdata))))
  changed = newdata[1:2, , drop = FALSE]
  changed[["twox"]][2L] = changed[["twox"]][2L] + 1
  expect_true(is.finite(predict(fit, changed)[[1L]]) && is.na(predict(fit, changed)[[2L]]))
  close(fit)
  # g enters as a numeric covariate here; all its variation is absorbed.
  fit = dbreg(y ~ x + s | g, conn = conn, table = table, weights = "w",
    strategy = "demean", retain = "prediction")
  models[[as.character(length(models) + 1L)]] = fit
  rows = newdata[1:3, , drop = FALSE]
  rows[["s"]][2L] = rows[["s"]][2L] + 1
  prediction = predict(fit, rows)
  expect_true(is.finite(prediction[[1L]]) && is.na(prediction[[2L]]) && is.finite(prediction[[3L]]))

  # The dictionary is independent of both outcome values and caller relations.
  expected = predict(fit, newdata)
  dbglm_validation_drop(conn, table)
  newdata[["y"]] = -999
  expect_equal(predict(fit, newdata), expected, tolerance = 1e-10)
  close(fit)
  expect_identical(dbglm_validation_tables(conn), sort(c(before, source)))
}
