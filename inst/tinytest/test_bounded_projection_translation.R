# Bounded FWL, SQL CG and iterative AP must use the same centered rank scale.
local({
  saved = options(dbreg.ap_solver = "auto", dbreg.ap_method = NULL)
  on.exit(options(saved), add = TRUE)
  set.seed(991)
  n = 600L
  d = data.frame(a = sample(1:12, n, TRUE), b = sample(1:13, n, TRUE),
    x = 1e14 + runif(n, -1, 1), z = rnorm(n), w = runif(n, .2, 2))
  d[["dx"]] = d[["x"]] - 1e14
  d[["absorbed"]] = 1e14 + d[["a"]] + d[["b"]]
  d[["y"]] = 2*d[["dx"]] + .3*d[["z"]] + d[["a"]] + d[["b"]] + rnorm(n, sd = .05)
  for (weighted in c(FALSE, TRUE)) {
    ref = lm(y ~ dx + z + factor(a) + factor(b), d,
      weights = if (weighted) w else NULL)
    for (strategy in c("auto", "demean", "compress")) {
      fit = dbreg(y ~ x + z + absorbed | a + b, data = d,
        weights = if (weighted) "w" else NULL, strategy = strategy)
      expect_identical(names(coef(fit)), c("x", "z"))
      expect_equal(unname(coef(fit)), unname(coef(ref)[c("dx", "z")]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[c("dx", "z"), c("dx", "z")]), tolerance = 1e-8)
    }
  }
})
