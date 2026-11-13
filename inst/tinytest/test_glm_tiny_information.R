# A nuisance group can have a finite identified mean even when its total prior
# weight makes its objective/ordinary standardized score negligible globally.
local({
  conn = DBI::dbConnect(duckdb::duckdb(config = list(temp_directory = "")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(g = c(rep(1L, 100), rep(2L, 3)),
    y = c(rep(c(0, 1), 50), 1, 1, 0), w = c(rep(1, 100), rep(1e-30, 3)))
  expected = c(rep(.5, 100), rep(2 / 3, 3))
  for (method in c("auto", "sql")) {
    f = dbglm(y ~ 1 | g, data = d, conn = conn, weights = ~w,
      family = binomial(), fixef_rm = "none", retain = "prediction",
      control = dbglm_control(projection_method = method))
    expect_true(f$converged && f$finite_mle)
    expect_equal(predict(f, d, type = "response"), expected, tolerance = 1e-9)
    expect_equal(predict(f, d), qlogis(expected), tolerance = 1e-8)
    expect_true(f$fe_score_max <= 1e-8)
    close(f)
    expect_identical(DBI::dbListTables(conn), character())
  }
  # The final equation's information guard uses base-group mass, so changing a
  # varying axis's units cannot create or erase its score certificate.
  working = data.frame(g = d$g, w = d$w / 4, score = d$w * (d$y - .5))
  DBI::dbWriteTable(conn, "working", working)
  score = dbreg:::dbglm_fe_group_score
  expect_equal(score(conn, "working", "g", nrow(d)), 2 / 3, tolerance = 1e-12)
  for (units in c(1e-10, 1, 1e10)) {
    ops = list(list(f = "g", v = sprintf("%.17g", units)))
    expect_equal(score(conn, "working", "g", nrow(d), ops), 2 / 3, tolerance = 1e-12)
  }
  DBI::dbRemoveTable(conn, "working")
  # An all-success negligible-weight group must not receive false finite-MLE
  # certification simply because its mass vanishes from the global objective.
  d$y[d$g == 2L] = 1
  expect_error(dbglm(y ~ 1 | g, data = d, conn = conn, weights = ~w,
    family = binomial(), fixef_rm = "none",
    control = dbglm_control(glm_max_iter = 8L)), "converge|boundary|finite|line|separation")
  expect_identical(DBI::dbListTables(conn), character())
})
