# Large running-variable locations must not enter geometry means, spline
# powers or covariance contrasts. References use the represented local data.
local({
  set.seed(571)
  n = 240L
  original = data.frame(x = 1e14 + runif(n, -1, 1))
  dx = original$x - 1e14
  original$z = 1 + sin(seq_len(n)*.37)
  original$g = rep(letters[1:4], length.out = n)
  original$cl = rep(seq_len(12L), length.out = n)
  original$y = 2*dx + .3*dx^2 + .4*original$z +
    .2*(match(original$g, letters[1:4])-1L) + rnorm(n, sd = .03)
  bin = pmin(4L, 1L + floor((dx-min(dx))/((max(dx)-min(dx))/4)))
  left = as.numeric(tapply(dx, bin, min))
  right = as.numeric(tapply(dx, bin, max))
  means = as.numeric(tapply(dx, bin, mean))
  knots = right[-4L]
  grid = unlist(Map(function(lo, hi) seq(lo, hi, length.out = 3), left, right))
  grid_bin = rep(1:4, each = 3)
  basis = function(x, bin, spec) {
    p = spec[1L]
    s = spec[2L]
    if (s == 0L) {
      X = outer(bin, 1:4, `==`)*1
      if (p > 0L) {
        u = (x-left[bin])/(right-left)[bin]
        for (power in seq_len(p)) X = cbind(X, outer(bin, 1:4, `==`)*u^power)
      }
    } else {
      X = matrix(1, length(x), 1)
      for (power in seq_len(p)) X = cbind(X, x^power)
      for (knot in knots) for (power in s:p) X = cbind(X, pmax(x-knot, 0)^power)
    }
    X
  }
  reference = function(spec, nuisance, kind) {
    X = basis(dx, bin, spec)
    Z = basis(means, 1:4, spec)
    line_Z = basis(grid, grid_bin, spec)
    if (nuisance) {
      extra = cbind(original$z, outer(original$g, letters[2:4], `==`)*1)
      X = cbind(X, extra)
      Z = cbind(Z, matrix(0, nrow(Z), ncol(extra)))
      line_Z = cbind(line_Z, matrix(0, nrow(line_Z), ncol(extra)))
    }
    qr_fit = lm.fit(X, original$y)
    if (qr_fit$rank != ncol(X)) stop("Reference design must have full column rank.")
    bread = solve(crossprod(X))
    e = qr_fit$residuals
    K = ncol(X)
    if (kind == "iid") V = bread*sum(e^2)/(n-K)
    else {
      score = X*e
      if (kind == "cluster") {
        score = rowsum(score, original$cl, reorder = FALSE)
        adjustment = 12/11*(n-1)/(n-K)
      } else adjustment = n/(n-K)
      V = adjustment*bread %*% crossprod(score) %*% bread
    }
    fit = drop(Z %*% qr_fit$coefficients)
    se = sqrt(diag(Z %*% V %*% t(Z)))
    critical = qt(.975, if (kind == "cluster") 11 else n-K)
    list(fit = fit, se = se, lwr = fit-critical*se, upr = fit+critical*se,
      line = drop(line_Z %*% qr_fit$coefficients))
  }
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  for (location in c(1e14, -1e14)) {
    d = original
    d$x = location + dx
    DBI::dbWriteTable(conn, "translated_bins", d, overwrite = TRUE)
    before = sort(DBI::dbListTables(conn))
    for (spec in list(c(0, 0), c(1, 0), c(1, 1), c(2, 2))) {
      for (nuisance in c(FALSE, TRUE)) for (kind in c("iid", "hc1", "cluster")) {
        result = suppressWarnings(dbbinsreg(if (nuisance) y ~ x + z | g else y ~ x,
          conn = conn, table = "translated_bins", nbins = 4, binspos = "es",
          points = spec, line = spec, linegrid = 3, plot = FALSE,
          vcov = if (kind == "cluster") ~cl else kind))
        expected = reference(spec, nuisance, kind)
        for (field in c("fit", "se", "lwr", "upr")) {
          expect_equal(result$points[[field]], expected[[field]], tolerance = 1e-7,
            info = paste(location, paste(spec, collapse = ","), nuisance, kind, field))
        }
        expect_equal(result$line$fit, expected$line, tolerance = 1e-7)
        expect_true(all(abs(result$points$x-(location+means)) <= abs(location)*.Machine$double.eps))
        expect_identical(result$bins$left, location+left)
        expect_identical(result$bins$right, location+right)
        if (spec[2L] > 0L) expect_identical(result$knots, location+knots)
        expect_identical(result$opt$x_var, "x")
        expect_identical(sort(DBI::dbListTables(conn)), before)
      }
    }
  }
  # Raw polynomial covariance already loses substantial precision at 1e4,
  # well before the most extreme representable-location example above.
  for (location in c(16, 1e4)) {
    d = original
    d$x = location + dx
    actual = suppressWarnings(dbbinsreg(y ~ x + z | g, data = d,
      nbins = 4, binspos = "es", points = c(2, 2), vcov = ~cl, plot = FALSE))
    expected = reference(c(2, 2), TRUE, "cluster")
    for (field in c("fit", "se", "lwr", "upr")) {
      expect_equal(actual$points[[field]], expected[[field]], tolerance = 1e-7)
    }
  }
})

local({
  # Ties, manual knots, separate point/line models and simulated bands retain
  # the same fitted quantities after translation; only reported x is restored.
  d = expand.grid(dx = seq(-1, 1, by = 1/8), rep = 1:8)
  d$x = 1e14 + d$dx
  d$y = 1 + 2*d$dx + .3*d$dx^2 + ifelse(d$rep %% 2, -.05, .05)
  for (position in list("qs", "es", 1e14 + c(-1, -.5, 0, .5, 1))) {
    reference_position = if (is.numeric(position)) position-1e14 else position
    set.seed(92)
    shifted = dbbinsreg(y ~ x, data = d, nbins = 4, binspos = position,
      points = c(1, 0), line = c(2, 2), linegrid = 5, cb = TRUE, nsims = 100, plot = FALSE)
    set.seed(92)
    centered = dbbinsreg(y ~ dx, data = d, nbins = 4, binspos = reference_position,
      points = c(1, 0), line = c(2, 2), linegrid = 5, cb = TRUE, nsims = 100, plot = FALSE)
    expect_equal(shifted$points[setdiff(names(shifted$points), "x")],
      centered$points[setdiff(names(centered$points), "x")], tolerance = 1e-8)
    expect_equal(shifted$line$fit, centered$line$fit, tolerance = 1e-8)
    expect_identical(shifted$bins$left, 1e14+centered$bins$left)
    expect_identical(shifted$bins$right, 1e14+centered$bins$right)
    expect_identical(shifted$knots$line, 1e14+centered$knots$line)
    expect_identical(shifted$opt$binspos, position)
  }
})

# An unused source column cannot shadow the private running-coordinate name.
local({
  set.seed(7441)
  d = data.frame(x = 1e14 + runif(240, -1, 1), dbbinsreg_running_x = rep(0, 240))
  d[["y"]] = 2*(d[["x"]]-1e14) + rnorm(nrow(d), sd = .1)
  fit = dbbinsreg(y ~ x, data = d, nbins = 4, binspos = "es", points = c(1, 1), plot = FALSE)
  d[["dbbinsreg_running_x"]] = NULL
  ref = dbbinsreg(y ~ x, data = d, nbins = 4, binspos = "es", points = c(1, 1), plot = FALSE)
  expect_equal(fit[["points"]], ref[["points"]], tolerance = 1e-10)
  expect_equal(fit[["opt"]][["x_origin"]], min(d[["x"]])/2+max(d[["x"]])/2)
})
