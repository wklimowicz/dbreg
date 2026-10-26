local({
  set.seed(1011)
  conn = DBI::dbConnect(duckdb::duckdb(), config = list(threads = "2"))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  resources = dbreg:::fit_resources()
  resources[["conn"]] = conn
  on.exit(dbreg:::release_fit_resources(resources), add = TRUE)
  manager = dbreg:::fit_table_manager(conn, resources)
  n = 3000L
  x = sample(-1:2, n, TRUE)
  a = sample(c(0.2, 0.5, 1, 1.5, 2.5), n, TRUE)
  offset = 0.1 * x
  for (family in list(stats::gaussian(), stats::poisson(), stats::binomial("logit"),
    stats::binomial("probit"), stats::binomial("cloglog"))) {
    eta = 0.3 * x + offset
    mu = family[["linkinv"]](eta)
    y = switch(family[["family"]], gaussian = mu + stats::rnorm(n),
      poisson = stats::rpois(n, mu), binomial = stats::rbinom(n, 1L, mu))
    source = data.frame(x = x)
    source[["__dbglm_key"]] = seq_len(n)
    source[["__dbglm_y"]] = y
    source[["__dbglm_a"]] = a
    source[["__dbglm_o"]] = offset
    source[["__dbglm_intercept"]] = 1
    DBI::dbWriteTable(conn, "cell_source", source, overwrite = TRUE)
    cells = dbreg:::dbglm_compress_base(conn, "cell_source",
      c("x", "__dbglm_o", "__dbglm_intercept"), family[["family"]],
      stats::weighted.mean(y, a), manager[["create"]], manager[["drop"]])
    X = cbind(1, x)
    derivative = family[["mu.eta"]](eta)
    variance = family[["variance"]](mu)
    w = a * derivative^2 / variance
    z = eta - offset + (y - mu) / derivative
    score = a * derivative * (y - mu) / variance
    terms = c(dbreg:::dbglm_moment_terms(c("__dbglm_intercept", "x"), "w"),
      "FSUM(w*z) AS b1", "FSUM(w*x*z) AS b2", "FSUM(dev) AS dev",
      "FSUM(score) AS s1", "FSUM(x*score) AS s2", "FSUM(score_sq) AS m11",
      "FSUM(x*score_sq) AS m12", "FSUM(x*x*score_sq) AS m22")
    sql = paste0(dbreg:::dbglm_working_sql(cells, c("__dbglm_intercept", "x"),
      c(0, 0.3), family[["link"]], cell = TRUE), " SELECT ", paste(terms, collapse = ","), " FROM working")
    row = dbreg:::backend_get_query(conn, sql)
    expect_equal(unname(dbreg:::dbglm_read_information(row, 2L)), unname(crossprod(X, w * X)), tolerance = 1e-11)
    expect_equal(as.numeric(row[1L, c("b1", "b2")]), as.numeric(crossprod(X, w * z)), tolerance = 1e-11)
    expect_equal(row[["dev"]], sum(family[["dev.resids"]](y, mu, a)), tolerance = 1e-11)
    expect_equal(as.numeric(row[1L, c("s1", "s2")]), as.numeric(crossprod(X, score)), tolerance = 1e-11)
    meat = crossprod(X, score^2 * X)
    expect_equal(as.numeric(row[1L, c("m11", "m12", "m22")]),
      c(meat[1L, 1L], meat[1L, 2L], meat[2L, 2L]), tolerance = 1e-11)
    if (family[["family"]] == "binomial") {
      totals = DBI::dbGetQuery(conn, paste("SELECT SUM(__dbglm_round_a) AS a,SUM(__dbglm_round_ay) AS ay FROM", cells))
      expect_equal(totals[["a"]], sum(round(a)))
      expect_equal(totals[["ay"]], sum(round(a) * y))
    }
    manager[["drop"]](cells)
  }
})
