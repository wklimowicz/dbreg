near = c(1, 1 + 2^-52, 1, 1 + 2^-51)
expect_identical(as.integer(dbreg:::exact_numeric_factor(near)), c(1L, 2L, 1L, 3L))
expect_identical(levels(dbreg:::exact_numeric_factor(c(10, 2, 1, 2))), c("1", "2", "10"))
expect_identical(levels(dbreg:::exact_integer_factor(c("10", "2", "-9007199254740992", "-9007199254740993"))),
  c("-9007199254740993", "-9007199254740992", "2", "10"))
large = bit64::as.integer64(c("9007199254740993", "-9007199254740993", "9007199254740992"))
expect_identical(as.integer(dbreg:::exact_numeric_factor(large)), c(3L, 1L, 2L))
