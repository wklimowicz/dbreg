# Safe-looking identifiers still need allocation against reserved targets.
expect_identical(dbreg:::safe_identifiers("x", taken = "x"), "x_1")
expect_identical(dbreg:::safe_identifiers("X", taken = "x"), "X_1")
expect_identical(dbreg:::safe_identifiers(c("x", "X", "x_1")), c("x", "X_2", "x_1"))
expect_identical(dbreg:::safe_identifiers(c("one column", "two column")), c("one_column", "two_column"))
expect_identical(dbreg:::safe_identifiers(c("a b", "a-b")), c("a_b", "a_b_1"))
expect_identical(dbreg:::safe_identifiers(c("a b", "a-b", "a_b")), c("a_b_1", "a_b_2", "a_b"))
expect_identical(dbreg:::safe_identifiers(character()), character())
expect_identical(dbreg:::safe_identifiers(c("A", "a"), fold_case = TRUE), c("a", "a_1"))

# A categorical dummy called grpb must not replace a raw regressor grpb.
local({
  set.seed(861)
  d = data.frame(y = rnorm(100), grp = rep(c("a", "b"), 50), grpb = rnorm(100))
  fit = dbreg(y ~ grp + grpb, data = d, strategy = "moments")
  ref = lm(y ~ grp + grpb, d)
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(ref)), tolerance = 1e-9)
  expect_equal(as.numeric(vcov(fit)), as.numeric(vcov(ref)), tolerance = 1e-9)
})
