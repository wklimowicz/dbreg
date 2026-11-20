# BIGINT aggregates must be decoded before entering ordinary double matrices.
# Exercise real DBI output rather than manufacturing an integer64 class.
local({
  modes = c("numeric", if (requireNamespace("bit64", quietly = TRUE)) "integer64")
  for (mode in modes) {
    conn = DBI::dbConnect(duckdb::duckdb(bigint = mode, config = list(threads = "2")))
    tryCatch({
      row = DBI::dbGetQuery(conn, paste(
        "SELECT CAST(10 AS BIGINT) AS sum_w, CAST(15 AS BIGINT) AS sum_w1,",
        "CAST(35 AS BIGINT) AS sum_w1_1, CAST(40 AS BIGINT) AS sum_w1_y,",
        "CAST(20 AS BIGINT) AS sum_wy, CAST(60 AS BIGINT) AS sum_wy_sq"))
      S = dbreg:::cross_moments(row, "x", "y", intercept = TRUE)
      expected = matrix(c(10, 15, 20, 15, 35, 40, 20, 40, 60), 3)
      expect_equal(as.numeric(S), as.numeric(expected))
      expect_equal(as.numeric(dbreg:::fit_ols(S, "y", c("(Intercept)", "x"))[["coef"]]), c(.8, .8))
      raw = DBI::dbGetQuery(conn, paste("SELECT CAST(3000000000 AS BIGINT) AS n,",
        "CAST(3000000001 AS BIGINT) AS dbreg_cell, CAST('9007199254740993' AS BIGINT) AS group_key"))
      decoded = dbreg:::decode_compressed_moments(raw)
      expect_identical(decoded[["n"]], 3e9)
      expect_identical(decoded[["dbreg_cell"]], 3000000001)
      expect_identical(decoded[["group_key"]], raw[["group_key"]])
    }, finally = DBI::dbDisconnect(conn, shutdown = TRUE))
  }
})

local({
  if (!requireNamespace("bit64", quietly = TRUE)) return(invisible(NULL))
  d = data.frame(x = rep(1:5, 60), cl = rep(1:12, each = 25))
  d$y = 2 + 0.5*d$x + sin(seq_len(nrow(d)))
  fit_mode = function(mode) {
    con = DBI::dbConnect(duckdb::duckdb(bigint = mode))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
    DBI::dbWriteTable(con, "d", d)
    lapply(c("moments", "compress", "auto"), function(strategy) {
      fit = dbreg(y ~ x, conn = con, table = "d", strategy = strategy, vcov = ~cl)
      list(beta = coef(fit), covariance = vcov(fit), df = fit$df_residual)
    })
  }
  expect_equal(fit_mode("integer64"), fit_mode("numeric"))
})
