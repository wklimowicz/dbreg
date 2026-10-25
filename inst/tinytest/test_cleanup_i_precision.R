local({
  set.seed(918)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (kind in c("numeric", "integer")) {
    keys = if (kind == "numeric") 1 + c(0, 1, 2) * .Machine[["double.eps"]] else
      c("9007199254740993", "9007199254740994", "9007199254740995")
    d = data.frame(id = 1:180, key = rep(keys, 60), index = factor(rep(1:3, 60)),
      x = rnorm(180), g = rep(1:9, each = 20))
    d[["y"]] = 0.7 * d[["x"]] + as.numeric(d[["index"]]) + sin(d[["g"]]) + rnorm(180, sd = 0.1)
    DBI::dbWriteTable(conn, "precision_source", d, overwrite = TRUE)
    DBI::dbExecute(conn, paste0("CREATE OR REPLACE TEMP VIEW precision_input AS SELECT id,x,y,g,",
      if (kind == "integer") "CAST(key AS BIGINT)" else "key", " AS key FROM precision_source"))
    for (fixed in c(FALSE, TRUE)) {
      fml = if (fixed) y ~ x + i(key) | g else y ~ x + i(key)
      reference = if (fixed) stats::lm(y ~ x + index + factor(g), d) else stats::lm(y ~ x + index, d)
      truth = as.numeric(stats::predict(reference, d))
      for (strategy in if (fixed) c("demean", "compress") else c("moments", "compress")) {
        for (retain in c("none", "all")) {
          fit = dbreg::dbreg(fml, conn = conn, table = "precision_input", strategy = strategy,
            retain = retain, row_id = "id")
          selected = c(if (!fixed) "(Intercept)", "x", "index2", "index3")
          expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)[selected]), tolerance = 1e-8)
          expect_equal(stats::predict(fit, d), truth, tolerance = 1e-8)
          expect_equal(length(fit[["i_specs"]][[1L]][["keys"]]), 3L)
          for (frame in if (!fixed || retain == "all") c(FALSE, TRUE) else logical()) {
            output = if (frame) stats::predict(fit, d, conn = conn, row_id = "id") else
              stats::predict(fit, conn = conn, table = "precision_input", row_id = "id")
            values = as.data.frame(output)
            expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])], truth, tolerance = 1e-8)
            close(output)
          }
          if (kind == "integer") {
            lossy = d
            lossy[["key"]] = as.numeric(lossy[["key"]])
            expect_error(stats::predict(fit, lossy), "must be exact")
            if (!fixed || retain == "all") expect_error(stats::predict(fit, lossy, conn = conn, row_id = "id"), "must be exact")
          }
          close(fit)
        }
      }
    }
  }
})

local({
  set.seed(921)
  d = data.frame(id = 1:120, x = rnorm(120), index = factor(rep(1:3, 40)), g = rep(1:6, each = 20))
  d[["key code"]] = rep(1 + c(0, 1, 2) * .Machine[["double.eps"]], 40)
  d[["y"]] = 0.5 * d[["x"]] + as.numeric(d[["index"]]) + sin(d[["g"]]) + rnorm(120, sd = 0.1)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (strategy in c("moments", "compress")) {
    fit = dbreg::dbreg(y ~ x + i(`key code`), data = d, strategy = strategy)
    missing = d
    missing[["key code"]][1:2] = c(NA_real_, Inf)
    truth = as.numeric(stats::predict(stats::lm(y ~ x + index, d), d))
    truth[1:2] = NA_real_
    expect_equal(stats::predict(fit, missing), truth, tolerance = 1e-8)
    output = stats::predict(fit, missing, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])], truth, tolerance = 1e-8)
    close(output)
    close(fit)
  }
  fit = dbreg::dbreg(y ~ x + i(`key code`) | g, data = d, retain = "prediction", strategy = "demean")
  mixed = d[1:3, ]
  mixed[["key code"]][1:2] = 99
  mixed[["g"]][1] = 99
  mixed[["x"]][2] = NA_real_
  truth = c(NA_real_, NA_real_, as.numeric(stats::predict(stats::lm(y ~ x + index + factor(g), d), d[3, ])))
  expect_equal(stats::predict(fit, mixed), truth, tolerance = 1e-8)
  output = stats::predict(fit, mixed, conn = conn, row_id = "id")
  values = as.data.frame(output)
  expect_equal(values[[output[["prediction_column"]]]][match(mixed[["id"]], values[["id"]])], truth, tolerance = 1e-8)
  close(output)
  close(fit)
  # A selector that removes every indicator still leaves a valid ordinary OLS design.
  fit = dbreg::dbreg(y ~ x + i(index, ref = c("1", "2", "3")), data = d, strategy = "moments")
  expect_equal(stats::predict(fit, d), as.numeric(stats::predict(stats::lm(y ~ x, d), d)), tolerance = 1e-8)
})

local({
  set.seed(930)
  d = data.frame(id = 1:180, y = rnorm(180), x = rnorm(180), key = rep(1:3, 60),
    f = rep(letters[1:3], 60), flag = rep(c(FALSE, TRUE), 90))
  cases = list(
    list(fml = y ~ x + key + i(key), reference = y ~ x + key + factor(key)),
    list(fml = y ~ x + f + i(f), reference = y ~ x + f),
    list(fml = y ~ x + flag + i(flag), reference = y ~ x + flag),
    list(fml = y ~ x + i(flag) + i(flag, x), reference = y ~ x * flag),
    list(fml = y ~ x + i(f, flag), reference = y ~ x + f:I(as.numeric(flag)))
  )
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (case in cases) for (strategy in c("moments", "compress")) {
    fit = dbreg::dbreg(case[["fml"]], data = d, strategy = strategy)
    reference = stats::lm(case[["reference"]], d)
    expect_equal(stats::predict(fit, d[2, , drop = FALSE]), as.numeric(stats::predict(reference, d[2, , drop = FALSE])), tolerance = 1e-8)
    output = stats::predict(fit, d, conn = conn, row_id = "id")
    values = as.data.frame(output)
    expect_equal(values[[output[["prediction_column"]]]][match(d[["id"]], values[["id"]])],
      as.numeric(stats::predict(reference, d)), tolerance = 1e-8)
    close(output)
    close(fit)
  }
})
