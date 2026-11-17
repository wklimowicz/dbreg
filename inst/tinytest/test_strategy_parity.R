# The same model gives the same answer whichever strategy or backend fits it,
# and that answer matches fixest.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")

set.seed(20261117)
n = 600
d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(20, n, TRUE),
  h = sample(8, n, TRUE), cl = sample(25, n, TRUE), w = runif(n, .5, 2))
d$y = .5*d$x - .2*d$z + rnorm(20)[d$g] + rnorm(8)[d$h] + rnorm(n)
full = fixest::ssc(K.fixef = "full")

## R2 is the total R2 for every strategy (demean used to report the within R2)
for (f in list(y ~ x + z | g, y ~ x + z | g + h)) for (wt in list(NULL, ~w)) {
  ref = fixest::feols(f, d, weights = wt, fixef.rm = "none")
  for (s in c("demean", "compress")) {
    fit = dbreg(f, data = d, strategy = s, weights = if (!is.null(wt)) "w")
    expect_equal(gof(fit)[["r2"]], unname(fixest::r2(ref, "r2")), tolerance = 1e-10,
      info = paste(deparse(f), s, !is.null(wt)))
  }
}
dm = data.frame(m = d$y, x = d$x, g = d$g)
expect_equal(gof(dbreg(m ~ x | g, data = dm, strategy = "demean"))[["r2"]],
  unname(fixest::r2(fixest::feols(m ~ x | g, dm, fixef.rm = "none"), "r2")), tolerance = 1e-10)
dw = data.frame(y = d$y, x = d$x, g = d$g, m = d$w)
expect_equal(gof(dbreg(y ~ x | g, data = dw, strategy = "demean", weights = "m"))[["r2"]],
  unname(fixest::r2(fixest::feols(y ~ x | g, dw, weights = ~m, fixef.rm = "none"), "r2")), tolerance = 1e-10)

## vcov() is a base matrix for every covariance type and strategy
for (s in c("moments", "compress")) for (kind in list("iid", "hc1", ~cl)) {
  fit = dbreg(y ~ x + z, data = d, strategy = s, vcov = kind, weights = "w")
  expect_true(is.matrix(vcov(fit)) && !isS4(vcov(fit)), info = paste(s, format(kind)))
  ref = fixest::feols(y ~ x + z, d, weights = ~w, vcov = kind, ssc = full)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-9,
    info = paste(s, format(kind)))
}
