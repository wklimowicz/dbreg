local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11169)
  core = expand.grid(g = 1:3, h = 1:3, x = -1:1, rep = 1:3)
  tail = data.frame(g = c(4, 4, 5, 5, 6), h = c(3, 4, 4, 5, 5), x = c(1, -1, 0, 1, -1), rep = 1)
  d = rbind(core, tail)
  n = nrow(d)
  d$w = runif(n, 0.3, 2)
  d$o = sin(d$rep) / 10
  d$cluster = rep(1:7, length.out = n)
  d$y = 1 + 0.4 * d$x + cos(d$g) + sin(d$h) + d$o + rnorm(n)
  d$x[1L] = NA
  d$w[2L] = 0
  d$id = seq_len(n)
  used = seq.int(3L, nrow(core))
  dense = lm(y ~ x + factor(g) + factor(h), data = d[used, ], weights = w, offset = o)
  for (strategy in c("demean", "compress")) {
    fit = dbreg(y ~ x | g + h, conn = conn, data = d, strategy = strategy,
      weights = ~w, offset = ~o, fixef_rm = "singletons")
    expect_equal(nobs(fit), length(used))
    expect_equal(fit$nobs_pruned, 5)
    expect_equal(fit$pruning_rounds, 5L)
    expect_identical(fit$fixef_rm, "singletons")
    expect_equal(coef(fit)[["x"]], coef(dense)[["x"]], tolerance = 1e-9)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(dense)["x", "x"]), tolerance = 1e-9)
    expect_identical(DBI::dbListTables(conn), character())
  }
  default = dbreg(y ~ x | g + h, conn = conn, data = d, strategy = "compress", weights = ~w, offset = ~o)
  expect_equal(nobs(default), n - 2L)
  expect_identical(default$fixef_rm, "none")
  expect_error(dbreg(y ~ x | g, conn = conn, data = d, fixef_rm = "all"), "arg")
  expect_error(dbreg(y ~ x | g, conn = conn, data = d, fixef_rm = "singletons", sql_only = TRUE), "requires execution")
  expect_error(dbreg(y ~ x | g, conn = conn, data = d, fixef_rm = "singletons", drop_missings = FALSE), "requires execution")
  expect_error(dbreg(y ~ x | id, conn = conn, data = d, fixef_rm = "singletons"), "No observations remain")
  expect_identical(DBI::dbListTables(conn), character())
})
