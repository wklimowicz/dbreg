local({
  set.seed(11011)
  d = data.frame(grp = rep(c("a", "b"), 60), grpb = stats::rnorm(120))
  d[["y"]] = 1 + 2 * (d[["grp"]] == "b") - 0.7 * d[["grpb"]] + stats::rnorm(120, sd = 0.2)
  rows = data.frame(id = c(17L, 4L, 90L, 2L, 31L, 8L),
    grp = c("a", "b", NA, "b", "a", "b"), grpb = c(2.5, -1.25, 0.3, NA, -2, 3.2))
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  for (factor_first in c(TRUE, FALSE)) {
    formula = if (factor_first) y ~ grp + grpb else y ~ grpb + grp
    design = function(data) {
      indicator = as.numeric(data[["grp"]] == "b")
      if (factor_first) cbind(1, indicator, data[["grpb"]]) else cbind(1, data[["grpb"]], indicator)
    }
    training_X = design(d)
    beta = as.numeric(qr.solve(training_X, d[["y"]]))
    residual = as.numeric(d[["y"]] - training_X %*% beta)
    reference_V = solve(crossprod(training_X)) * sum(residual^2) / (nrow(d) - ncol(training_X))
    valid = !is.na(rows[["grp"]]) & !is.na(rows[["grpb"]])
    expected = rep(NA_real_, nrow(rows))
    expected[valid] = as.numeric(design(rows[valid, ]) %*% beta)
    for (strategy in c("moments", "compress")) {
      fit = dbreg::dbreg(formula, conn = conn, data = d, strategy = strategy)
      expect_equal(as.numeric(stats::coef(fit)), beta, tolerance = 1e-10)
      expect_equal(as.numeric(stats::vcov(fit)), as.numeric(reference_V), tolerance = 1e-10)
      expect_equal(stats::df.residual(fit), nrow(d) - ncol(training_X))
      expect_equal(attr(stats::vcov(fit), "df.K"), ncol(training_X))
      expect_equal(stats::nobs(fit), nrow(d))
      expect_true(anyDuplicated(names(stats::coef(fit))) > 0L)
      expect_equal(stats::predict(fit, rows), expected, tolerance = 1e-10)
      # Both interval kinds use the independent full positional covariance.
      new_X = design(rows[valid, , drop = FALSE])
      mean_variance = rowSums((new_X %*% reference_V) * new_X)
      sigma_squared = sum(residual^2) / (nrow(d) - ncol(training_X))
      critical = stats::qt(0.975, nrow(d) - ncol(training_X))
      for (interval in c("confidence", "prediction")) {
        width = critical * sqrt(mean_variance + if (interval == "prediction") sigma_squared else 0)
        reference_interval = data.frame(fit = expected, lwr = NA_real_, upr = NA_real_)
        reference_interval[valid, "lwr"] = expected[valid] - width
        reference_interval[valid, "upr"] = expected[valid] + width
        expect_equal(stats::predict(fit, rows, interval = interval), reference_interval, tolerance = 1e-10)
      }
      changed_type = rows[1:2, ]
      changed_type[["grpb"]] = factor(c("100", "200"))
      # The grp dummy still has the public label grpb. Its presence must not
      # disguise a numeric-source type change in the other grpb coordinate.
      expect_error(stats::predict(fit, changed_type), "numeric term|rebuild the design")
      extra = rows
      extra[["dbreg_internal_x_0001__"]] = seq_len(nrow(rows)) * 100
      extra[["irrelevant"]] = NA_real_
      expect_equal(stats::predict(fit, extra), expected, tolerance = 1e-10)
      shuffled = rev(seq_len(nrow(rows)))
      expect_equal(stats::predict(fit, extra[shuffled, ]), expected[shuffled], tolerance = 1e-10)
      for (table_source in c(FALSE, TRUE)) {
        if (table_source) {
          DBI::dbWriteTable(conn, "native_column_rows", extra[shuffled, ], overwrite = TRUE)
          out = stats::predict(fit, conn = conn, table = "native_column_rows", row_id = "id")
        } else out = stats::predict(fit, extra[shuffled, ], conn = conn, row_id = "id")
        actual = as.data.frame(out)
        expect_equal(sort(actual[["id"]]), sort(rows[["id"]]))
        expect_equal(actual[[out[["prediction_column"]]]][match(rows[["id"]], actual[["id"]])],
          expected, tolerance = 1e-10)
        close(out)
      }
      close(fit)
    }
  }
})
