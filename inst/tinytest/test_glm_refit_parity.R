# A prepared refit equals a fresh dbglm fit with the same weights.
local({
  set.seed(20261117)
  n = 600
  d = data.frame(id = seq_len(n), x = rnorm(n), z = rnorm(n), g = sample(15, n, TRUE), w = runif(n, .5, 2))
  d$y = rpois(n, exp(.3*d$x - .2*d$z + rnorm(15, 0, .3)[d$g]))
  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  DBI::dbWriteTable(conn, "d", d)
  prep = dbglm_prepare(y ~ x + z | g, conn = conn, table = "d", row_id = "id", family = poisson(),
    weights = "rw", etastart = "eta0", vcov = "iid")
  for (k in 1:3) {
    nw = if (k == 1) d$w else runif(n, .2, 3)
    DBI::dbWriteTable(conn, "state", data.frame(id = d$id, rw = nw, eta0 = 0), overwrite = TRUE)
    rf = dbglm_refit(prep, "state", etastart = FALSE)
    DBI::dbWriteTable(conn, "fresh", transform(d, rw = nw), overwrite = TRUE)
    fresh = dbglm(y ~ x + z | g, conn = conn, table = "fresh", family = poisson(), weights = "rw", vcov = "iid")
    expect_equal(coef(rf), coef(fresh), tolerance = 1e-9, info = k)
    expect_equal(as.matrix(vcov(rf)), as.matrix(vcov(fresh)), tolerance = 1e-8, info = k)
  }
  close(prep)
})
