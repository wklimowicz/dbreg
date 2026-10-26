local({
  set.seed(43001)
  d = data.frame(x = stats::runif(160))
  d[["y"]] = 2 + 3 * d[["x"]] + 0.7 * d[["x"]]^2 + 0.1 * sin(17 * d[["x"]])
  breaks = seq(0, 1, length.out = 5)
  bins = as.integer(cut(d[["x"]], breaks, include.lowest = TRUE, right = FALSE))
  left = as.numeric(tapply(d[["x"]], bins, min))
  right = as.numeric(tapply(d[["x"]], bins, max))
  knots = right[1:3]
  design = function(x, bin, degree, smooth) {
    if (!smooth) {
      b = factor(bin, levels = seq_len(4))
      intercepts = stats::model.matrix(~b)
      indicators = stats::model.matrix(~b - 1)
      u = (x - left[bin]) / (right[bin] - left[bin])
      powers = lapply(seq_len(degree), function(power) indicators * u^power)
      do.call(cbind, c(list(intercepts), powers))
    } else {
      columns = c(list(rep(1, length(x))), lapply(seq_len(degree), function(power) x^power))
      for (knot in knots) for (power in smooth:degree) columns = c(columns, list(pmax(x - knot, 0)^power))
      do.call(cbind, columns)
    }
  }
  for (spec in list(c(1, 0), c(2, 0), c(1, 1), c(2, 2))) {
    X = design(d[["x"]], bins, spec[1], spec[2])
    reference = stats::lm.fit(X, d[["y"]])
    actual = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
      points = FALSE, line = spec, ci = FALSE, sample_fit = FALSE, plot = FALSE)
    grid = actual[["line"]]
    expected = as.numeric(design(grid[["x"]], grid[["bin"]], spec[1], spec[2]) %*% reference[["coefficients"]])
    expect_equal(grid[["fit"]], expected, tolerance = 1e-8)
    standalone = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
      points = spec, line = spec, ci = FALSE, sample_fit = FALSE, plot = FALSE)
    expect_equal(grid, standalone[["line"]], tolerance = 1e-10)
    expect_identical(actual[["points"]], NULL)
    expect_identical(actual[["opt"]][["line"]], as.integer(spec))
    expect_identical(actual[["model"]][["vcov_type"]], "iid")
  }
  implicit = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = FALSE, line = TRUE, ci = FALSE, plot = FALSE)
  explicit = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = FALSE, line = c(0, 0), ci = FALSE, plot = FALSE)
  expect_equal(implicit[["line"]], explicit[["line"]], tolerance = 1e-12)
  neither = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = FALSE, line = FALSE, ci = FALSE, plot = FALSE)
  expect_identical(neither[["points"]], NULL)
  expect_identical(neither[["line"]], NULL)
  expect_true(!is.null(neither[["model"]][["coeftable"]]))
})
