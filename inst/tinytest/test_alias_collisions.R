expect_equal(dbreg:::safe_identifiers("x", taken = "x"), "x_1")
expect_equal(dbreg:::safe_identifiers(c("x", "X", "x_1")), c("x", "X_2", "x_1"))
expect_equal(dbreg:::safe_identifiers(c("one column", "two column")), c("one_column", "two_column"))

# A syntactically valid dummy alias must still avoid an existing regressor.
set.seed(861)
dat = data.frame(y = rnorm(100), grp = rep(c("a", "b"), 50), grpb = rnorm(100))
fit = dbreg(y ~ grp + grpb, data = dat, strategy = "moments")
ref = lm(y ~ grp + grpb, dat)
expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
