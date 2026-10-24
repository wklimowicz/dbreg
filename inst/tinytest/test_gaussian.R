source(system.file("tinytest/helper_glm_validation.R", package = "dbreg"), local = TRUE)
source(system.file("tinytest/helper_gaussian.R", package = "dbreg"), local = TRUE)
local({
  conn = dbglm_validation_connection()
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))
  dbreg_test_gaussian(conn)
  set.seed(569)
  d = data.frame(x = rnorm(144), t = rnorm(144), g = rep(1:8, each = 18), h = rep(1:6, 24), o = rnorm(144))
  d[["y"]] = .4*d[["x"]]+d[["g"]]*d[["t"]]/8+d[["o"]]+rnorm(144)
  for (fml in list(y ~ x | g[t]+h, y ~ 0 | g, y ~ 1)) {
    fit = dbglm(fml, data = d, family = gaussian(), offset = "o", retain = "all", fixef_rm = "none",
      control = dbglm_control(projection_method = "sql"))
    ref = fixest::feglm(fml, d, family = "gaussian", offset = ~o, fixef.rm = "none", notes = FALSE)
    expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-7)
    expect_equal(predict(fit, d), predict(ref, d), tolerance = 1e-7)
    expect_equal(predict(fit), predict(ref), tolerance = 1e-7)
    expect_equal(residuals(fit), residuals(ref), tolerance = 1e-7)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-7)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-7)
    close(fit)
  }
})
