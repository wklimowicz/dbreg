# Translating a predictor cannot erase representable within-FE variation.
# Center the reference explicitly so its QR never sees the large location.
local({
  saved = options(dbreg.ap_solver = "iterative", dbreg.ap_tol = 1e-9)
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
    for (method in c("join", "window")) {
      options(dbreg.ap_method = method)
      fit = dbreg(y ~ x + z + absorbed | a + b, data = d,
        weights = if (weighted) "w" else NULL, strategy = "demean")
      expect_identical(names(coef(fit)), c("x", "z"))
      expect_identical(fit[["collin.var"]], "absorbed")
      expect_equal(unname(coef(fit)), unname(coef(ref)[c("dx", "z")]), tolerance = 1e-8)
      expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)[c("dx", "z"), c("dx", "z")]), tolerance = 1e-8)
    }
  }
})
