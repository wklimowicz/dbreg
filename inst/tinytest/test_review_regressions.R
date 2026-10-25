# Public value invariants and failure cleanup from the whole-package review.
local({
  set.seed(440)
  n = 1000L
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = rep(1:20, length.out = n))
  d[["y"]] = 2+.5*d[["x"]]-.3*d[["z"]]+rnorm(n)
  d[["w"]] = runif(n, .3, 2)
  for (strategy in c("moments", "compress")) {
    message("Review invariants: ", strategy)
    for (column in c("x", "y")) {
      shifted = d
      shifted[[column]] = shifted[[column]] + if (column == "x") 1e6 else 1e8
      actual = dbreg(y ~ x+z, data = shifted, strategy = strategy)
      reference = lm(y ~ x+z, shifted)
      expect_equal(unname(coef(actual)), unname(coef(reference)), tolerance = 1e-6)
      expect_equal(as.numeric(vcov(actual)), as.numeric(vcov(reference)), tolerance = 1e-6)
      expect_equal(attr(vcov(actual), "rss"), sum(residuals(reference)^2), tolerance = 1e-7)
      for (kind in list("hc1", ~g)) {
        baseline = dbreg(y ~ x+z, data = d, strategy = strategy, vcov = kind)
        robust = dbreg(y ~ x+z, data = shifted, strategy = strategy, vcov = kind)
        A = diag(3)
        if (column == "x") A[1,2] = -1e6
        expected = A %*% vcov(baseline) %*% t(A)
        expect_true(all(abs(vcov(robust)-expected) <= 1e-8+1e-6*abs(expected)))
      }
    }
    for (kind in list("iid", "hc1", ~g)) {
      baseline = dbreg(y ~ x+z, data = d, weights = "w", vcov = kind, strategy = strategy)
      for (scale in c(1e-200, 1e200)) {
        scaled = d
        scaled[["w"]] = scaled[["w"]]*scale
        fit = dbreg(y ~ x+z, data = scaled, weights = "w", vcov = kind, strategy = strategy)
        expect_equal(coef(fit), coef(baseline), tolerance = 1e-8)
        expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(baseline)), tolerance = 1e-8)
      }
    }
  }
  # Integer products must be promoted before multiplication, including default selection.
  integer = data.frame(x = 50000L+seq_len(n), y = seq_len(n)+rnorm(n))
  expect_equal(unname(coef(dbreg(y ~ x, data = integer))), unname(coef(lm(y ~ x, integer))), tolerance = 1e-7)
  expect_error(dbreg(y ~ x, data = data.frame(x = 1:2, y = c(2.1, 4.2)), strategy = "moments"), pattern = "degrees of freedom")
  one = d
  one[["g"]] = 1
  expect_error(dbreg(y ~ x, data = one, vcov = ~g, strategy = "moments"), pattern = "at least two clusters")

  conn = DBI::dbConnect(duckdb::duckdb())
  message("Review invariants: missing values and inference")
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  # NaNs are missing values in both R and the native SQL sample.
  for (column in c("y", "x", "w", "g")) {
    missing = d
    missing[[column]][1:6] = NaN
    ref = lm(y ~ x+z, missing, weights = w, subset = !is.na(g))
    actual = dbreg(y ~ x+z, data = missing, weights = "w", vcov = ~g, strategy = "moments")
    expect_equal(nobs(actual), nobs(ref))
    expect_equal(unname(coef(actual)), unname(coef(ref)), tolerance = 1e-8)
    DBI::dbWriteTable(conn, "review_sample", missing, overwrite = TRUE)
    # Preserve a real SQL NaN rather than depending on upload conversion.
    DBI::dbExecute(conn, paste0("UPDATE review_sample SET ", column, "=CAST('NaN' AS DOUBLE) WHERE id<=6"))
    sql = dbreg(y ~ x+z, conn = conn, table = "review_sample", weights = "w", vcov = ~g, strategy = "moments")
    expect_equal(nobs(sql), nobs(ref))
    expect_equal(coef(sql), coef(actual), tolerance = 1e-8)
  }
  few = d
  few[["g"]] = rep(1:3, length.out = n)
  fit = dbreg(y ~ x+z, data = few, conn = conn, vcov = ~g, strategy = "moments", retain = "all", row_id = "id")
  on.exit(close(fit), add = TRUE)
  recomputed = summary(fit, vcov = ~g)
  expect_equal(as.numeric(fit[["coeftable"]]), as.numeric(recomputed[["coeftable"]]), tolerance = 1e-8)
  expect_equal(as.numeric(confint(fit)), as.numeric(confint(recomputed)), tolerance = 1e-8)
  expect_equal(attr(vcov(fit), "df.t"), 2)
  for (strategy in c("moments", "compress")) {
    translated = few
    translated[["x"]] = translated[["x"]]+1e6
    translated[["y"]] = translated[["y"]]+1e8
    retained = dbreg(y ~ x+z, data = translated, conn = conn, vcov = ~g, strategy = strategy, retain = "all", row_id = "id")
    V = vcov(retained)
    updated = vcov(retained, vcov = ~g)
    expect_true(all(abs(V-updated) <= 1e-8+1e-6*abs(V)))
    close(retained)
  }

  # Observation-sized cluster dictionaries have the HC1 covariance under CR1.
  unique = d
  unique[["g"]] = unique[["id"]]
  robust = dbreg(y ~ x+z, data = unique, vcov = "hc1", strategy = "moments")
  clustered = dbreg(y ~ x+z, data = unique, vcov = ~g, strategy = "moments")
  expect_equal(as.numeric(vcov(clustered)), as.numeric(vcov(robust)), tolerance = 1e-9)

  # SQL binscatter agrees across actual database backends, with negative x and controls.
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    message("Review invariants: binscatter backends")
    sqlite = DBI::dbConnect(RSQLite::SQLite(), ":memory:")
    on.exit(DBI::dbDisconnect(sqlite), add = TRUE)
    DBI::dbWriteTable(sqlite, "review_bins", d)
    DBI::dbWriteTable(conn, "review_bins", d)
    for (degree in 0:2) {
      a = dbbinsreg(y ~ x+z, conn = conn, table = "review_bins", points = c(degree, 0), nbins = 5, vcov = "hc1", plot = FALSE)
      b = dbbinsreg(y ~ x+z, conn = sqlite, table = "review_bins", points = c(degree, 0), nbins = 5, vcov = "hc1", plot = FALSE)
      expect_equal(a[["points"]], b[["points"]], tolerance = 1e-7)
    }
  }
  for (strategy in c("auto", "compress")) {
    base = dbbinsreg(y ~ x, data = d, nbins = 5, strategy = strategy, vcov = "hc1", plot = FALSE)
    translated = d
    translated[["y"]] = translated[["y"]]+1e8
    shifted = dbbinsreg(y ~ x, data = translated, nbins = 5, strategy = strategy, vcov = "hc1", plot = FALSE)
    expect_equal(shifted[["points"]][["fit"]]-1e8, base[["points"]][["fit"]], tolerance = 1e-6)
    expect_equal(shifted[["points"]][["se"]], base[["points"]][["se"]], tolerance = 1e-6)
    if (strategy == "compress") expect_equal(base[["model"]][["strategy"]], "compress")
  }
  before = sort(DBI::dbListTables(conn))
  invalid = data.frame(x = factor(rep(letters[1:3], 4000)), y = rnorm(12000))
  expect_error(dbbinsreg(y ~ x, data = invalid, conn = conn, points = c(1,1), nbins = 3, randcut = .5, sample_fit = TRUE, plot = FALSE))
  expect_identical(sort(DBI::dbListTables(conn)), before)
})
