# A small absolute remainder of an FE-only column must not turn into a
# seemingly identified regressor when the ordinary regressor is rescaled.
local({
  saved = options(dbreg.ap_tol = 1e-7, dbreg.ap_max_iter = 200L,
    dbreg.ap_method = getOption("dbreg.ap_method"))
  on.exit(options(saved), add = TRUE)
  set.seed(431)
  d = data.frame(x = rnorm(400), g = sample(1:10, 400, TRUE), h = sample(1:7, 400, TRUE), w = runif(400, .2, 2))
  d[["absorbed"]] = sin(d[["g"]]) + cos(d[["h"]])
  d[["y"]] = .4*d[["x"]] + d[["absorbed"]] + rnorm(400)
  conn = DBI::dbConnect(duckdb::duckdb(config = list(threads = "2")))
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  for (weighted in c(FALSE, TRUE)) {
    ref = fixest::feols(y ~ x | g + h, d, weights = if (weighted) ~w else NULL,
      vcov = "hc1", fixef.rm = "none", fixef.tol = 1e-10,
      ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
    for (method in c("window", "join")) {
      options(dbreg.ap_method = method)
      for (scale in c(1, 1e-6, 1e-12)) {
        scaled = d
        scaled[["x"]] = scaled[["x"]]*scale
        fit = dbreg(y ~ x + absorbed | g + h, data = scaled, conn = conn,
          strategy = "demean", weights = if (weighted) "w" else NULL, vcov = "hc1")
        expect_identical(names(coef(fit)), "x")
        expect_equal(as.numeric(coef(fit))*scale, as.numeric(coef(ref)), tolerance = 1e-7)
        expect_equal(as.numeric(vcov(fit))*scale^2, as.numeric(vcov(ref)), tolerance = 1e-7)
        expect_identical(DBI::dbListTables(conn), character())
      }
    }
  }
  # Tight iteration caps fail honestly and release all AP state.
  options(dbreg.ap_max_iter = 1L)
  for (method in c("join", "window")) {
    options(dbreg.ap_method = method)
    expect_error(dbreg(y ~ x | g + h, data = d, conn = conn,
      strategy = "demean", weights = "w"), "did not converge")
    expect_true(DBI::dbIsValid(conn))
    expect_identical(DBI::dbListTables(conn), character())
  }
})
