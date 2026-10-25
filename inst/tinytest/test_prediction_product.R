local({
  multiply = dbreg:::prediction_matrix_product
  set.seed(931)
  x = matrix(rnorm(500),100,5)
  beta = c(.1,2,-3,0,1)
  expect_identical(multiply(x,beta),as.vector(x%*%beta))
  expect_equal(multiply(Matrix::Matrix(x,sparse=TRUE),beta),as.vector(x%*%beta),tolerance=1e-14)
  extreme = rbind(c(1,1e308,1e308),c(1,1e308,0),c(1,-1e308,0),c(1,0,0))
  value = multiply(extreme,c(1,2,-2))
  expect_identical(is.nan(value),c(TRUE,FALSE,FALSE,FALSE))
  expect_identical(is.infinite(value),c(FALSE,TRUE,TRUE,FALSE))
  expect_true(is.nan(multiply(extreme[1,,drop=FALSE],c(1,2,-2))))
  expect_identical(value[[4L]],1)
  # Partial-sum overflow follows SQL's staged arithmetic, without extended
  # accumulators or contraction depending on BLAS and matrix batch dimensions.
  expect_identical(multiply(matrix(1,1,3),c(1e308,1e308,-1e308)),Inf)
  expect_identical(multiply(matrix(1,1,3),c(-1e308,1e308,1e308)),1e308)
  expect_identical(multiply(matrix(numeric(),4,0),numeric()),rep(0,4))
})
