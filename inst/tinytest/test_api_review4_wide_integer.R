# Native learned FE keys can be wider than signed BIGINT. SQL upload and
# comparisons must retain the same exact domain as portable R prediction.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(654)
  domains = list(
    decimal = list(type = "DECIMAL(38,0)", keys = c(
      "-99999999999999999999999999999999999999", "-100000000000000000000",
      "-1", "0", "1", "100000000000000000000", "99999999999999999999999999999999999999")),
    huge = list(type = "HUGEINT", keys = c(
      "-170141183460469231731687303715884105728", "-9223372036854775809",
      "-1", "0", "1", "9223372036854775808", "170141183460469231731687303715884105727")),
    bigint = list(type = "BIGINT", keys = c("-9223372036854775807", "-9007199254740993",
      "-1", "0", "1", "9007199254740993", "9223372036854775807")))
  for (domain in domains) {
    keys = domain[["keys"]]
    d = data.frame(id = seq_len(25L*length(keys)), x = rnorm(25L*length(keys)),
      group_number = rep(seq_along(keys), each = 25L))
    d[["y"]] = .7*d[["x"]] + sin(d[["group_number"]]) + rnorm(nrow(d), sd = .1)
    DBI::dbWriteTable(conn, "review4_integer_source", d, overwrite = TRUE)
    values = paste0("(", seq_along(keys), ",CAST('", keys, "' AS ", domain[["type"]], "))")
    DBI::dbExecute(conn, paste0("CREATE OR REPLACE TEMP TABLE review4_integer_keys AS SELECT s.*,k.g FROM review4_integer_source s JOIN (VALUES ",
      paste(values, collapse = ","), ") k(group_number,g) ON s.group_number=k.group_number"))
    fit = dbreg(y ~ x | g, conn = conn, table = "review4_integer_keys", row_id = "id", retain = "all")
    expected = as.numeric(fitted(lm(y ~ x + factor(group_number), d)))
    text = data.frame(id = d[["id"]], x = d[["x"]], g = keys[d[["group_number"]]])
    expect_equal(predict(fit, text), expected, tolerance = 1e-9)
    sql = predict(fit, conn = conn, table = "review4_integer_keys", row_id = "id")
    collected = as.data.frame(sql)
    expect_equal(collected[[sql[["prediction_column"]]]], expected, tolerance = 1e-9)
    close(sql)
    # Keep the raw output text, including accepted signs and leading zeroes.
    text[["g"]] = ifelse(startsWith(text[["g"]], "-"),
      paste0("-000", substring(text[["g"]], 2L)), paste0("+000", text[["g"]]))
    sql = predict(fit, text, row_id = c("id", "g"))
    collected = as.data.frame(sql)
    expect_identical(collected[["g"]], text[["g"]])
    expect_equal(collected[[sql[["prediction_column"]]]], expected, tolerance = 1e-9)
    close(sql)
    if (domain[["type"]] == "BIGINT") {
      text[["g"]] = bit64::as.integer64(keys[d[["group_number"]]])
      expect_equal(predict(fit, text), expected, tolerance = 1e-9)
      sql = predict(fit, text, row_id = c("id", "g"))
      collected = as.data.frame(sql)
      expect_identical(as.character(collected[["g"]]), as.character(text[["g"]]))
      expect_equal(collected[[sql[["prediction_column"]]]], expected, tolerance = 1e-9)
      close(sql)
    }
    bad = data.frame(id = 1:6, x = 0, g = c("1.0", "1e0", " 1", "bad", "42",
      "170141183460469231731687303715884105728"))
    expect_true(all(is.na(predict(fit, bad))))
    sql = predict(fit, bad, row_id = "id")
    expect_true(all(is.na(as.data.frame(sql)[[sql[["prediction_column"]]]])))
    close(sql)
    # Native fractional DECIMAL inputs must not round into integer FE levels.
    DBI::dbExecute(conn, "CREATE OR REPLACE TEMP TABLE review4_fractional_keys AS SELECT * FROM (VALUES (1,0.0,CAST(0.1 AS DECIMAL(38,1))),(2,0.0,CAST(0.5 AS DECIMAL(38,1))),(3,0.0,CAST(-0.5 AS DECIMAL(38,1))),(4,0.0,CAST(1 AS DECIMAL(38,1)))) t(id,x,g)")
    sql = predict(fit, conn = conn, table = "review4_fractional_keys", row_id = "id")
    collected = as.data.frame(sql)[[sql[["prediction_column"]]]]
    expect_true(all(is.na(collected[1:3])))
    expect_equal(collected[4L], as.numeric(predict(fit, data.frame(x = 0, g = "1"))), tolerance = 1e-9)
    close(sql)
    close(fit)
  }
})
