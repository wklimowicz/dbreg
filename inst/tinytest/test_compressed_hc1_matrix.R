# The same model gives the same answer whichever strategy or backend fits it,
# and that answer matches fixest.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")

set.seed(20261117)
n = 600
d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(20, n, TRUE),
  h = sample(8, n, TRUE), cl = sample(25, n, TRUE), w = runif(n, .5, 2))
d$y = .5*d$x - .2*d$z + rnorm(20)[d$g] + rnorm(8)[d$h] + rnorm(n)
full = fixest::ssc(K.fixef = "full")

## vcov() is a base matrix for every covariance type and strategy
for (s in c("moments", "compress")) for (kind in list("iid", "hc1", ~cl)) {
  fit = dbreg(y ~ x + z, data = d, strategy = s, vcov = kind, weights = "w")
  expect_true(is.matrix(vcov(fit)) && !isS4(vcov(fit)), info = paste(s, format(kind)))
  ref = fixest::feols(y ~ x + z, d, weights = ~w, vcov = kind, ssc = full)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-9,
    info = paste(s, format(kind)))
}
