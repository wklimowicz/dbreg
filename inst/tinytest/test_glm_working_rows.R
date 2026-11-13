local({
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  d = data.frame(x = seq(-2, 2, length.out = 21L),
    `__dbglm_y` = rep(c(0, 1, 1), 7L),
    `__dbglm_a` = seq(.25, 2.75, length.out = 21L),
    `__dbglm_o` = seq(.1, -.2, length.out = 21L), check.names = FALSE)
  table = "glm_working_rows"
  DBI::dbWriteTable(conn, table, d)
  families = list(gaussian(), poisson(), binomial("logit"), binomial("probit"),
    binomial("cloglog"))
  for (family in families) {
    link = family[["link"]]
    rows = dbreg:::dbglm_working_rows_sql(table, "x", .7, link, conn = conn)
    got = DBI::dbGetQuery(conn, paste0(rows, "SELECT * FROM working ORDER BY x"))
    eta = d[["__dbglm_o"]] + .7 * d[["x"]]
    mu = family[["linkinv"]](eta)
    derivative = family[["mu.eta"]](eta)
    variance = family[["variance"]](mu)
    y = d[["__dbglm_y"]]
    weight = d[["__dbglm_a"]]
    expect_equal(got[["mu"]], mu, tolerance = 1e-12)
    expect_equal(got[["w"]], weight * derivative^2 / variance, tolerance = 1e-12)
    expect_equal(got[["score"]], weight * derivative * (y - mu) / variance,
      tolerance = 1e-12)
    expect_equal(got[["dev"]], family[["dev.resids"]](y, mu, weight),
      tolerance = 1e-12)
    expect_equal(got[["z"]], eta - d[["__dbglm_o"]] + (y - mu) / derivative,
      tolerance = 1e-12)
    if (family[["family"]] == "binomial") {
      observed = dbreg:::dbglm_binomial_loglik(conn, paste0("(", rows,
        "SELECT * FROM working) likelihood"))
      expect_equal(observed, sum(round(weight) * dbinom(y, 1, mu, log = TRUE)),
        tolerance = 1e-12)
    }
  }
})
