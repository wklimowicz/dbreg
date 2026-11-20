# Covariance extraction follows the coefficients exposed by coef(). The
# explicit-FE matrix is still required for full-design prediction intervals.
local({
  set.seed(810)
  d = data.frame(x = rnorm(240), f = factor(rep(1:8, 30)), cl = rep(1:12, each = 20))
  d[["y"]] = 2*d[["x"]] + as.numeric(d[["f"]]) + rnorm(nrow(d))
  for (kind in list("iid", "hc1", ~cl)) {
    fit = dbreg(y ~ x | f, data = d, strategy = "compress", vcov = kind)
    full = fit[["vcov"]]
    visible = names(coef(fit))
    expect_identical(rownames(vcov(fit)), visible)
    expect_identical(colnames(vcov(fit)), visible)
    expect_equal(as.numeric(vcov(fit)), as.numeric(full[visible, visible, drop = FALSE]))
    expect_equal(vcov(fit, fe = TRUE), full)
    metadata = c("type", "n_clusters", "n_params", "rss", "tss")
    for (name in metadata) expect_identical(attr(vcov(fit), name), attr(full, name))
    oracle = fixest::feols(y ~ x | f, data = d, vcov = kind, fixef.rm = "none",
      ssc = fixest::ssc(K.fixef = "full", K.exact = TRUE))
    expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(oracle)), tolerance = 1e-9)
  }
  fit = dbreg(y ~ x | f, data = d, strategy = "compress")
  ref = lm(y ~ x + f, d)
  for (interval in c("confidence", "prediction")) {
    expect_equal(as.numeric(as.matrix(predict(fit, d, interval = interval))),
      as.numeric(predict(ref, d, interval = interval)), tolerance = 1e-9)
  }
  fit = dbreg(y ~ x * f, data = d, strategy = "compress")
  expect_identical(rownames(vcov(fit)), names(coef(fit)))
  expect_identical(colnames(vcov(fit)), names(coef(fit)))
})

# Display labels may repeat. Without FE filtering, extract by position rather
# than matching names, which would silently repeat the first parameter.
local({
  V = matrix(c(2, .3, .3, 5), 2, dimnames = list(c("a", "b"), c("a", "b")))
  attr(V, "type") = "iid"
  attr(V, "rss") = 7
  ct = cbind(estimate = c(.4, .8))
  rownames(ct) = c("same", "same")
  fit = structure(list(vcov = V, coeftable = ct, fe = NULL), class = "dbreg")
  expect_equal(as.numeric(vcov(fit)), as.numeric(V))
  expect_identical(rownames(vcov(fit)), c("same", "same"))
  expect_identical(colnames(vcov(fit)), c("same", "same"))
  expect_identical(attr(vcov(fit), "type"), "iid")
  expect_identical(attr(vcov(fit), "rss"), 7)
})
