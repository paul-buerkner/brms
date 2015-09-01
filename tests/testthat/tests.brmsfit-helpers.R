test_that("Test that first_greater returns expected results", {
  set.seed(12345)
  A <- cbind(1:10, 11:20, 21:30)
  x <- sample(1:25, size = 10, replace = TRUE)
  expect_equal(first_greater(A, x), c(3,3,3,3,2,1,2,2,2,3))
  expect_equal(first_greater(A, x, i = 2), c(3,3,3,3,2,2,2,2,2,3))
})

test_that("Test that array2list performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(matrix(1:9,3,3), matrix(10:18,3,3), matrix(19:27,3,3))
  expect_equal(array2list(A), B)
})

test_that("Test that list2array performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(array(1:9, dim = c(3,3,1)), array(10:18, dim = c(3,3,1)), array(19:27, dim = c(3,3,1)))
  expect_equal(list2array(B), A)
})