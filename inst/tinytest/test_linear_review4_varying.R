# Shifted structural variation must survive compressed varying-slope rank
# checks. Compare to a full, stable weighted dummy/trend row design.
local({
  set.seed(76376)
  n = 180L
  d = data.frame(g = rep(1:6, each = 30), h = rep(1:3, length.out = n),
    q = rnorm(n), trend = 1e12 + rnorm(n), w = runif(n, .4, 1.3))
  d[["x"]] = 1e14 + runif(n, -1, 1)
  d[["dx"]] = d[["x"]] - 1e14
  d[["dt"]] = d[["trend"]] - 1e12
  d[["joint"]] = interaction(d[["g"]], d[["h"]], drop = TRUE)
  d[["y"]] = 2*d[["dx"]] + .4*d[["q"]] + sin(d[["g"]]) + .2*d[["g"]]*d[["dt"]] + rnorm(n, sd = .05)
  saved = options(dbreg.ap_solver = "auto")
  on.exit(options(saved))
  for (joint in c(FALSE, TRUE)) {
    ref = lm(if (joint) y ~ dx + q + joint + joint:dt else y ~ dx + q + factor(g) + factor(g):dt,
      d, weights = w)
    fml = if (joint) y ~ x + q | (g^h)[trend] else y ~ x + q | g[trend]
    fit = dbreg(fml, data = d, weights = "w", strategy = "compress")
    expect_identical(names(coef(fit)), c("x", "q"))
    expect_equal(unname(coef(fit)), unname(coef(ref)[c("dx", "q")]), tolerance = 1e-8)
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[c("dx", "q"), c("dx", "q")]), tolerance = 1e-8)
    expect_equal(as.numeric(predict(fit, d)), as.numeric(fitted(ref)), tolerance = 1e-8)
  }
  # Genuine null axes and absorbed real controls still disappear.
  d[["absorbed"]] = 1e14 + d[["g"]]
  d[["constant"]] = 1
  fit = dbreg(y ~ x + q + absorbed | g[trend, constant], data = d,
    weights = "w", strategy = "compress")
  expect_identical(names(coef(fit)), c("x", "q"))
  expect_true("absorbed" %in% fit[["collin.var"]])
})
