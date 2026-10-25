# Disconnected graphs have additional dummy-matrix rank deficiencies. Solver
# representation must not change the public full-FE SSC counting convention.
local({
  set.seed(5727)
  n = 360L
  component = rep(0:2, each = n/3)
  d = data.frame(x = rnorm(n), z = rnorm(n),
    g = component*4L + rep(1:4, length.out = n),
    h = component*2L + sample.int(2L, n, TRUE), cl = rep(1:12, length.out = n))
  d[["y"]] = .4*d[["x"]]-.2*d[["z"]]+sin(d[["g"]])+rnorm(n)
  for (kind in c("iid", "hc1", "cluster")) {
    vc = if (kind == "cluster") ~cl else kind
    reference = fixest::feols(y ~ x + z | g + h, d,
      vcov = if (kind == "hc1") "hetero" else vc,
      fixef.rm = "none", ssc = fixest::ssc(K.fixef = "full"))
    for (strategy in c("compress", "demean")) {
      fit = dbreg(y ~ x + z | g + h, data = d, strategy = strategy, vcov = vc)
      expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(reference)), tolerance = 1e-8)
      close(fit)
    }
  }
})
