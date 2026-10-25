# Typed original-row SQL scores replace the former client-side cell join.
local({
  set.seed(719)
  d = expand.grid(a = c("a\rb", "a"), b = c("c", "b\rc"), cl = 1:4, replicate = 1:20,
    stringsAsFactors = FALSE)
  d[["x"]] = rnorm(nrow(d))
  d[["y"]] = .4*d[["x"]]+(d[["a"]] == "a")+rnorm(nrow(d))
  a = dbreg(y ~ x+a+b, data = d, strategy = "compress", vcov = ~cl)
  b = dbreg(y ~ x+a+b, data = d, strategy = "moments", vcov = ~cl)
  expect_equal(coef(a), coef(b), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(a)), as.numeric(vcov(b)), tolerance = 1e-8)
  d[["f"]] = rep(c(1, 1+2^-52), length.out = nrow(d))
  d[["cl"]] = rep(rep(c(1, 1+2^-52), each = 3), length.out = nrow(d))
  a = dbreg(y ~ x | f, data = d, strategy = "compress", vcov = ~cl)
  b = dbreg(y ~ x | f, data = d, strategy = "demean", vcov = ~cl)
  expect_equal(coef(a), coef(b), tolerance = 1e-8)
  expect_equal(as.numeric(vcov(a)), as.numeric(vcov(b)), tolerance = 1e-8)
  expect_equal(attr(vcov(a), "n_clusters"), 2)
})
