# User-specified geometry is the oracle, including gaps around interior knots.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  d = data.frame(x = c(seq(.05, .25, length.out = 30), seq(.75, .95, length.out = 30)))
  d[["y"]] = 2 + d[["x"]] + 8 * pmax(d[["x"]] - .5, 0) + .03 * sin(1:60)
  reference = stats::lm(y ~ x + I(pmax(x - .5, 0)), data = d)
  for (position in list(c(0, .5, 1), "es")) {
    fit = dbbinsreg(y ~ x, conn = conn, data = d, binspos = position, nbins = 2,
      points = c(1, 1), line = c(1, 1), linegrid = 3, vcov = "iid", plot = FALSE)
    expect_equal(as.numeric(fit[["knots"]]), .5, tolerance = 1e-14)
    expected = predict(reference, newdata = data.frame(x = fit[["points"]][["x"]]), se.fit = TRUE)
    expect_equal(fit[["points"]][["fit"]], unname(expected[["fit"]]), tolerance = 1e-9)
    expect_equal(fit[["points"]][["se"]], unname(expected[["se.fit"]]), tolerance = 1e-9)
    expect_equal(fit[["line"]][["fit"]], unname(predict(reference,
      newdata = data.frame(x = fit[["line"]][["x"]]))), tolerance = 1e-9)
  }
  # An unusable-outcome row cannot change the equal-width partition.
  d = data.frame(x = c(100, 0:9), y = c(NA, (0:9)^2))
  fit = dbbinsreg(y ~ x, data = d, conn = conn, nbins = 2, binspos = "es", plot = FALSE)
  expect_equal(fit[["points"]][["x"]], c(2, 7))
  expect_equal(fit[["points"]][["fit"]], c(6, 51))
  expect_equal(fit[["opt"]][["N"]], 10)

  # Ordinary sample quantiles retain their own exact endpoints in both routes.
  set.seed(2026)
  d = data.frame(x = runif(120))
  d[["y"]] = 1 + 2 * d[["x"]] + sin(1:120)
  breaks = unname(quantile(d[["x"]], seq(0, 1, length.out = 5), type = 2))
  bin = findInterval(d[["x"]], breaks, rightmost.closed = TRUE)
  means = as.numeric(tapply(d[["x"]], bin, mean))
  for (smooth in 0:1) {
    fit = dbbinsreg(y ~ x, data = d, conn = conn, binspos = breaks,
      points = c(smooth, smooth), plot = FALSE, vcov = "iid")
    expect_equal(nobs(fit[["model"]]), nrow(d))
    expect_equal(fit[["opt"]][["N"]], nrow(d))
    expect_equal(fit[["points"]][["x"]], means, tolerance = 1e-14)
    if (smooth == 0L) {
      expect_equal(fit[["points"]][["fit"]], as.numeric(tapply(d[["y"]], bin, mean)), tolerance = 1e-10)
    } else {
      knots = breaks[2:4]
      X = cbind(1, d[["x"]], vapply(knots, function(k) pmax(d[["x"]] - k, 0), numeric(nrow(d))))
      Z = cbind(1, means, vapply(knots, function(k) pmax(means - k, 0), numeric(length(means))))
      expected = stats::lm.fit(X, d[["y"]])
      expect_equal(fit[["knots"]], knots)
      expect_equal(fit[["points"]][["fit"]], drop(Z %*% expected[["coefficients"]]), tolerance = 1e-9)
    }
  }
  breaks = c(0, .50000000000000044, 1)
  d = data.frame(x = rep(c(.1, .3, .5, .50000000000000022, .7, .9), each = 5))
  d[["y"]] = rep(1:6, each = 5) + rep(c(-.02, -.01, 0, .01, .02), 6)
  fit = dbbinsreg(y ~ x, data = d, conn = conn, binspos = breaks, plot = FALSE)
  expect_equal(fit[["points"]][["fit"]], c(2.5, 5.5), tolerance = 1e-12)
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})

# All fitted quantities agree with an explicitly filtered estimation sample.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  initial = sort(DBI::dbListTables(conn))
  i = 1:80
  original = data.frame(x = (i - .5)/80, z = sin(i * .71),
    g = rep(letters[1:4], 20), cl = rep(1:10, 8))
  original[["y"]] = 1 + original[["x"]] + 2 * pmax(original[["x"]] - .5, 0) +
    .7 * original[["z"]] + sin(i * 1.17)/10
  for (variable in c("z", "g", "cl")) {
    d = original
    d[[variable]][c(1:20, 41:50)] = NA
    complete = d[stats::complete.cases(d), ]
    for (position in list(c(0, .5, 1), "qs", "es")) {
      for (spec in list(c(0, 0), c(1, 1))) {
        args = list(fml = y ~ x + z | g, conn = conn, binspos = position,
          nbins = 2, points = spec, line = spec, linegrid = 3, vcov = ~cl, plot = FALSE)
        fit = do.call(dbbinsreg, c(args, list(data = d)))
        reference = do.call(dbbinsreg, c(args, list(data = complete)))
        expect_equal(fit[["points"]], reference[["points"]], tolerance = 1e-9)
        expect_equal(fit[["line"]], reference[["line"]], tolerance = 1e-9)
        expect_equal(fit[["bins"]], reference[["bins"]], tolerance = 1e-12)
        expect_equal(fit[["knots"]], reference[["knots"]], tolerance = 1e-12)
        expect_equal(fit[["opt"]][["N"]], nrow(complete))
        expect_equal(nobs(fit[["model"]]), nrow(complete))
      }
    }
    separated = dbbinsreg(y ~ x + z | g, conn = conn, data = d, nbins = 2,
      points = c(0, 0), line = c(1, 1), linegrid = 3, vcov = ~cl, plot = FALSE)
    reference = dbbinsreg(y ~ x + z | g, conn = conn, data = complete, nbins = 2,
      points = c(0, 0), line = c(1, 1), linegrid = 3, vcov = ~cl, plot = FALSE)
    expect_equal(separated[["points"]], reference[["points"]], tolerance = 1e-9)
    expect_equal(separated[["line"]], reference[["line"]], tolerance = 1e-9)
    expect_equal(separated[["knots"]], reference[["knots"]], tolerance = 1e-12)
    expect_true(all(vapply(separated[["model"]], nobs, numeric(1)) == nrow(complete)))
  }
  expect_identical(sort(DBI::dbListTables(conn)), initial)
})
