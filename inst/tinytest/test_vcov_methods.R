set.seed(810)
dat = data.frame(x = rnorm(240), f = factor(rep(1:8, 30)), cl = rep(1:12, each = 20))
dat[["y"]] = 2 * dat[["x"]] + as.numeric(dat[["f"]]) + rnorm(nrow(dat))
for (kind in list("iid", "hc1", ~cl)) {
  fit = dbreg(y ~ x | f, data = dat, strategy = "compress", vcov = kind)
  full = fit[["vcov"]]
  expect_equal(rownames(vcov(fit)), names(coef(fit)))
  expect_equal(colnames(vcov(fit)), names(coef(fit)))
  expect_equal(as.matrix(vcov(fit)), as.matrix(full["x", "x", drop = FALSE]), check.attributes = FALSE)
  expect_equal(vcov(fit, fe = TRUE), full)
  expect_equal(attr(vcov(fit), "rss"), attr(full, "rss"))
  expect_equal(attr(vcov(fit), "type"), attr(full, "type"))
}
fit = dbreg(y ~ x | f, data = dat, strategy = "compress")
ref = lm(y ~ x + f, dat)
expect_equal(
  as.matrix(predict(fit, dat, interval = "confidence")),
  predict(ref, dat, interval = "confidence"), check.attributes = FALSE, tolerance = 1e-8
)
interaction_fit = dbreg(y ~ x * f, data = dat, strategy = "compress")
expect_equal(rownames(vcov(interaction_fit)), names(coef(interaction_fit)))
expect_equal(colnames(vcov(interaction_fit)), names(coef(interaction_fit)))

# R permits a dummy label to coincide with a different raw regressor name.
# With no FE to hide, covariance extraction must preserve parameter positions.
dat[["grpb"]] = rnorm(nrow(dat))
dat[["grp"]] = rep(c("a", "b"), length.out = nrow(dat))
fit = dbreg(y ~ grp + grpb, data = dat, strategy = "moments")
ref = lm(y ~ grp + grpb, dat)
expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-8)
