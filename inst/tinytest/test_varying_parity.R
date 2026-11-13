# Varying slopes match fixest's g[z] for every strategy.
if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
local({
  set.seed(20261117)
  n = 800
  d = data.frame(x = rnorm(n), z = rnorm(n), g = sample(20, n, TRUE), h = sample(8, n, TRUE),
    cl = sample(25, n, TRUE), w = runif(n, .5, 2))
  d$y = .5*d$x + rnorm(20)[d$g] + rnorm(20, .3)[d$g]*d$z + rnorm(8)[d$h] + rnorm(n)
  full = fixest::ssc(K.fixef = "full")
  for (f in list(y ~ x | g[z], y ~ x | g[z] + h)) for (v in list("iid", "hc1", ~cl)) for (s in c("auto", "compress")) {
    info = paste(deparse(f), format(v), s)
    ref = fixest::feols(f, d, weights = ~w, vcov = v, ssc = full, fixef.tol = 1e-11, fixef.rm = "none")
    fit = dbreg(f, data = d, weights = "w", vcov = v, strategy = s)
    expect_equal(coef(fit)[["x"]], coef(ref)[["x"]], tolerance = 1e-8, info = info)
    expect_equal(sqrt(diag(as.matrix(vcov(fit))))[["x"]], fixest::se(ref)[["x"]], tolerance = 1e-7, info = info)
  }
})

## Varying slopes beyond the dense budget use coordinate SQL CG, which needs more than the
## 100 alternating-projection cycles on a low-mobility worker-firm panel.
set.seed(7)
nw = 300; nf = 120; T = 6
wk = rep(seq_len(nw), each = T); t = rep(seq_len(T), nw)
home = sample(nf, nw, TRUE); move = runif(nw) < 0.3; dest = sample(nf, nw, TRUE); mt = sample(2:T, nw, TRUE)
fm = ifelse(move[wk] & t >= mt[wk], dest[wk], home[wk])
n = length(wk); x = rnorm(n); z = rnorm(n)
y = 0.5*x - 0.2*z + rnorm(nw)[wk] + rnorm(nw, 0, .05)[wk]*t + rnorm(nf)[fm] + rnorm(n)
d = data.frame(y, x, z, wk, fm, t, cl = fm %% 40)
f = y ~ x + z | wk[t] + fm
ref = fixest::feols(f, d, vcov = ~cl, ssc = fixest::ssc(K.fixef = "full"), fixef.tol = 1e-11,
  fixef.iter = 1e5, fixef.rm = "none")
old = options(dbreg.fe_dense_max = 0)
fit = dbreg(f, data = d, vcov = ~cl)
options(old)
expect_equal(coef(fit)[c("x", "z")], coef(ref)[c("x", "z")], tolerance = 1e-8)
expect_equal(sqrt(diag(vcov(fit)))[c("x", "z")], fixest::se(ref)[c("x", "z")], tolerance = 1e-6)
