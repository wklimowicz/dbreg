# Reproduce the remaining translation sensitivity in raw-moment rank detection.
pkgload::load_all()
set.seed(7)
base = data.frame(x = rnorm(1000))
base[["y"]] = 2 + 3 * base[["x"]] + rnorm(1000)
results = lapply(c(0, 100, 1000, 10000, 1000000), function(shift) {
  d = base
  d[["x"]] = d[["x"]] + shift
  got = dbreg(y ~ x, data = d, strategy = "moments")
  ref = lm(y ~ x, d)
  data.frame(shift = shift, retained_x = "x" %in% names(coef(got)),
    dbreg_slope = if ("x" %in% names(coef(got))) coef(got)[["x"]] else NA_real_,
    lm_slope = coef(ref)[["x"]], rss = attr(vcov(got), "rss"))
})
print(do.call(rbind, results), row.names = FALSE)
