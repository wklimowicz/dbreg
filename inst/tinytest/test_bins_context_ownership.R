# Saved metadata must not retain the caller or the fitting frame.
local({
  make_result = function(bytes) {
    sentinel = raw(bytes)
    set.seed(42001)
    d = data.frame(x = stats::runif(80), cluster = rep(seq_len(8), 10))
    d[["y"]] = 1 + d[["x"]] + stats::rnorm(80)
    request = ~cluster
    original_environment = environment(request)
    original_call = as.call(request)
    result = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = seq(0, 1, length.out = 5),
      points = c(0, 0), vcov = request, plot = FALSE)
    expect_identical(environment(request), original_environment)
    expect_identical(as.call(request), original_call)
    expect_identical(environment(result[["opt"]][["formula"]]), baseenv())
    expect_identical(environment(result[["opt"]][["vcov"]]), baseenv())
    expect_identical(environment(result[["model"]][["fml"]]), baseenv())
    result
  }
  small = make_result(10000L)
  large = make_result(1000000L)
  expect_equal(small[["points"]], large[["points"]], tolerance = 1e-12)
  # serialize() follows non-global formula environments; object.size() does not.
  expect_true(abs(length(serialize(large, NULL)) - length(serialize(small, NULL))) < 10000L)
  restored = unserialize(serialize(large, NULL))
  expect_equal(restored[["points"]], large[["points"]], tolerance = 1e-12)
  expect_identical(environment(restored[["opt"]][["formula"]]), baseenv())
})

# Composition retains each branch's numerical result and effective point policy.
local({
  set.seed(42002)
  d = data.frame(x = stats::runif(160))
  d[["y"]] = 2 + 3 * d[["x"]] + stats::rnorm(160)
  breaks = seq(0, 1, length.out = 5)
  points = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(0, 0), ci = TRUE, plot = FALSE)
  line = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(1, 1), line = c(1, 1), ci = FALSE, sample_fit = FALSE, plot = FALSE)
  combined = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(0, 0), line = c(1, 1), ci = TRUE, sample_fit = FALSE, plot = FALSE)
  expect_equal(combined[["points"]], points[["points"]], tolerance = 1e-10)
  expect_equal(combined[["line"]], line[["line"]], tolerance = 1e-10)
  expect_identical(combined[["opt"]][["ci"]], TRUE)
  expect_identical(combined[["opt"]][["cb"]], FALSE)
  expect_identical(combined[["model"]][["line"]][["vcov_type"]], "iid")
  expect_identical(combined[["opt"]][["points"]], c(0L, 0L))
  expect_identical(combined[["opt"]][["line"]], c(1L, 1L))

  set.seed(42003)
  bands = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(0, 0), cb = TRUE, nsims = 30, plot = FALSE)
  set.seed(42003)
  combined_bands = dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(0, 0), line = c(1, 1), cb = TRUE, nsims = 30, sample_fit = FALSE, plot = FALSE)
  expect_equal(combined_bands[["points"]], bands[["points"]], tolerance = 1e-10)
  expect_identical(combined_bands[["opt"]][["cb"]], TRUE)
  expect_warning(dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(1, 1), line = c(0, 0), cb = TRUE, sample_fit = FALSE, plot = FALSE), "not yet supported")
  unsupported = suppressWarnings(dbreg::dbbinsreg(y ~ x, data = d, nbins = 4, binspos = breaks,
    points = c(1, 1), line = c(0, 0), cb = TRUE, sample_fit = FALSE, plot = FALSE))
  expect_identical(unsupported[["opt"]][["cb"]], FALSE)
  expect_identical(unsupported[["opt"]][["ci"]], TRUE)
})

# Outer cleanup remains the only owner of borrowed caller connection/source.
local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(42004)
  d = data.frame(x = stats::runif(80), y = stats::rnorm(80))
  DBI::dbWriteTable(conn, "caller_source", d)
  sparse = data.frame(x = rep(c(.1, .35, .6, .85), 20), y = stats::rnorm(80))
  DBI::dbWriteTable(conn, "caller_sparse", sparse)
  tables = sort(DBI::dbListTables(conn))
  result = dbreg::dbbinsreg(y ~ x, conn = conn, table = "caller_source", nbins = 4,
    binspos = seq(0, 1, length.out = 5), points = c(0, 0), line = c(1, 1), sample_fit = FALSE, plot = FALSE)
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), tables)
  expect_equal(DBI::dbReadTable(conn, "caller_source"), d)
  expect_error(dbreg::dbbinsreg(y ~ x, conn = conn, table = "caller_sparse", nbins = 4,
    binspos = seq(0, 1, length.out = 5), points = c(0, 0), line = c(3, 0), plot = FALSE), "enough distinct")
  expect_true(DBI::dbIsValid(conn))
  expect_identical(sort(DBI::dbListTables(conn)), tables)
  expect_equal(DBI::dbReadTable(conn, "caller_sparse"), sparse)
})
