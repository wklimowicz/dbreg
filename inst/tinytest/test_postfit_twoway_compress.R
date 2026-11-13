# Post-fit two-way covariance is repaired on the reported coefficients, as in fixest.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(9); n = 3000
  d = data.frame(id = 1:n, x = rbinom(n, 1, .4), z = sample(0:4, n, TRUE), firm = sample(40, n, TRUE),
    state = sample(25, n, TRUE), year = sample(12, n, TRUE), f = sample(letters[1:8], n, TRUE))
  d$y = d$x - .5 * d$z + d$firm / 30 + rnorm(n)
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  ref = fixest::feols(y ~ x + z | firm, d, fixef.tol = 1e-11)
  se_ref = fixest::se(summary(ref, cluster = ~state + year, ssc = fixest::ssc(K.fixef = "full")))
  for (s in c("compress", "demean")) {
    fit = dbreg(y ~ x + z | firm, conn = conn, table = "d", strategy = s, retain = "inference", row_id = "id")
    expect_silent(V <- vcov(fit, vcov = ~state + year))
    expect_equal(as.numeric(sqrt(diag(V))), as.numeric(se_ref), tolerance = 1e-8, info = s)
    sm = summary(fit, vcov = ~state + year)
    expect_equal(as.numeric(sm[["coeftable"]][c("x", "z"), "std.error"]), as.numeric(se_ref), tolerance = 1e-8, info = s)
    close(fit)
  }
  ref2 = fixest::feols(y ~ x + z + factor(f), d)
  se_ref2 = suppressWarnings(fixest::se(summary(ref2, cluster = ~state + year)))
  fit2 = dbreg(y ~ x + z + factor(f), conn = conn, table = "d", strategy = "compress", retain = "inference", row_id = "id")
  se2 = suppressWarnings(sqrt(diag(vcov(fit2, vcov = ~state + year))))
  expect_equal(as.numeric(se2), as.numeric(se_ref2), tolerance = 1e-8)
  # Intervals use the same repaired covariance as the reported coefficients.
  nd = d[1:100, ]
  pref = suppressWarnings(predict(summary(ref2, cluster = ~state + year), newdata = nd, interval = "confidence"))
  sm2 = suppressWarnings(summary(fit2, cluster = ~state + year))
  p2 = as.matrix(suppressWarnings(predict(sm2, newdata = nd, interval = "confidence")))
  expect_equal(as.numeric(p2[, 3] - p2[, 1]), pref$ci_high - pref$fit, tolerance = 1e-8)
  close(fit2)
  ref3 = suppressWarnings(summary(fixest::feols(y ~ x + z + factor(firm), d), cluster = ~state + year,
    ssc = fixest::ssc(K.fixef = "full")))
  pref3 = predict(ref3, newdata = nd, interval = "confidence")
  fit3 = dbreg(y ~ x + z | firm, conn = conn, table = "d", strategy = "compress", retain = "inference", row_id = "id")
  p3 = as.matrix(suppressWarnings(predict(suppressWarnings(summary(fit3, cluster = ~state + year)), newdata = nd, interval = "confidence")))
  expect_equal(as.numeric(p3[, 3] - p3[, 1]), pref3$ci_high - pref3$fit, tolerance = 1e-8)
  close(fit3)
})
