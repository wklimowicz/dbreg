# Entire FE spans, not just individual FE terms, must remain absorbed.
local({
  saved = options(dbreg.ap_solver = "iterative", dbreg.ap_tol = 1e-7, dbreg.ap_max_iter = 200L, dbreg.ap_method = NULL)
  on.exit(options(saved), add = TRUE)
  set.seed(431)
  d = data.frame(x = rnorm(400), g = sample(1:10, 400, TRUE), h = sample(1:7, 400, TRUE), w = runif(400, .2, 2))
  d[["absorbed"]] = sin(d[["g"]])+cos(d[["h"]])
  d[["y"]] = .4*d[["x"]]+d[["absorbed"]]+rnorm(400)
  ref = fixest::feols(y ~ x | g+h, d, weights = ~w, vcov = "hc1", fixef.rm = "none",
    fixef.tol = 1e-10, ssc = fixest::ssc(K.fixef = "full"), notes = FALSE)
  for (method in c("window", "join")) {
    options(dbreg.ap_method = method)
    for (scale in c(1, 1e-6)) {
      scaled = d
      scaled[["x"]] = scaled[["x"]]*scale
      fit = dbreg(y ~ x+absorbed | g+h, data = scaled, strategy = "demean", weights = "w", vcov = "hc1")
      expect_identical(names(coef(fit)), "x")
      expect_equal(coef(fit)*scale, coef(ref), tolerance = 1e-7)
      expect_equal(as.numeric(vcov(fit))*scale^2, as.numeric(vcov(ref)), tolerance = 1e-7)
    }
  }
})
