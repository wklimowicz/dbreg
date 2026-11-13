local({
  conn = DBI::dbConnect(duckdb::duckdb(bigint = "integer64"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  set.seed(11167)
  d = data.frame(g = rep(1:5, each = 20), x = rep(-1:1, length.out = 100),
    w = runif(100, .1, 2), batch = rep(1:2, 50))
  d$y = rpois(100, exp(.2 + .3*d$x + cos(d$g)/4))
  d$y[d$g == 5L] = 0
  d$x[4L] = NA
  d$w[5L] = 0
  d[["id key"]] = bit64::as.integer64("9007199254740992") + rep(1:50, each = 2)
  DBI::dbWriteTable(conn, "source", d, temporary = TRUE)
  before = DBI::dbReadTable(conn, "source")
  selected = with(d, g != 5L & !is.na(x) & w > 0 & batch == 1L)
  expected = d[selected, c("id key", "batch")]
  dense = glm(y ~ x + factor(g), data = d[selected, ], family = poisson(), weights = w)
  for (strategy in c("moments", "compress")) {
    for (retain in c("none", "all")) {
      model = dbglm(y ~ x | g, conn = conn, table = "source", family = poisson(),
        weights = ~w, subset = ~batch == 1L, row_id = c("id key", "batch"),
        sample_table = "export", strategy = strategy, retain = retain)
      expect_identical(model$sample_table, "export")
      expect_equal(nobs(model), nrow(expected))
      expect_equal(unname(coef(model)), unname(coef(dense)["x"]), tolerance = 1e-8)
      close(model)
      actual = DBI::dbGetQuery(conn, 'SELECT CAST("id key" AS VARCHAR) AS "id key",batch FROM export ORDER BY "id key",batch')
      expect_identical(names(actual), names(expected))
      expect_identical(as.character(actual[["id key"]]), as.character(expected[["id key"]]))
      expect_equal(actual$batch, expected$batch)
      DBI::dbExecute(conn, "DROP TABLE export")
      expect_identical(DBI::dbListTables(conn), "source")
    }
  }
  expect_error(dbglm(y ~ x, conn = conn, table = "source", sample_table = "generated"), "row_id")
  expect_error(dbglm(y ~ x, conn = conn, table = "source", row_id = "id key", sample_table = "export"), "unique and nonmissing")
  expect_error(dbglm(y ~ x, data = d, row_id = c("id key", "batch"), sample_table = "export"), "caller-owned")
  expect_error(dbglm(y ~ x, conn = conn, table = "source", row_id = c("id key", "batch"), sample_table = ""), "sample_table")
  expect_error(dbglm(y ~ x | g, conn = conn, table = "source", family = poisson(), weights = ~w,
    row_id = c("id key", "batch"), sample_table = "failed", control = dbglm_control(glm_max_iter = 1L)), "converge")
  expect_false(DBI::dbExistsTable(conn, "failed"))
  expect_error(dbglm(y ~ x, conn = conn, table = "source", row_id = "batch"), "unique and nonmissing")
  expect_identical(DBI::dbListTables(conn), "source")
  expect_identical(DBI::dbReadTable(conn, "source"), before)
})
