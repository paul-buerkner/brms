test_that("student distribution works correctly", {
  expect_equal(integrate(dstudent, -100, 100, df = 15, mu = 10, sigma = 5)$value, 1)
  expect_equal(dstudent(1, df = 10, mu = 0, sigma = 5), dt(1/5, df = 10)/5)
  expect_equal(pstudent(2, df = 20, mu = 2, sigma = 0.4), pt(0, df = 20))
  expect_equal(qstudent(0.7, df = 5, mu = 2, sigma = 3), 2 + 3*qt(0.7, df = 5))
  expect_equal(length(rstudent(10, df = 10, mu = rnorm(10), sigma = 1:10)), 10)
})

test_that("multivariate normal and student distributions work correctly", {
  library(mvtnorm)
  mu <- rnorm(3)
  Sigma <- cov(matrix(rnorm(300), ncol = 3))
  expect_equal(dmulti_normal(1:3, mu = mu, Sigma = Sigma),
               dmvnorm(1:3, mean = mu, sigma = Sigma, log = TRUE))
  expect_equal(dmulti_student(1:3, mu = mu, Sigma = Sigma, df = 10),
               dmvt(1:3, df = 10, delta = mu, sigma = Sigma, log = TRUE))
  expect_equal(dim(rmulti_normal(7, mu = mu, Sigma = Sigma)), c(7, 3))
  expect_equal(dim(rmulti_student(7, mu = mu, Sigma = Sigma, df = 10)), 
               c(7, 3))
  # test errors
  expect_error(dmulti_normal(1:3, mu = rnorm(2), Sigma = Sigma, check = TRUE),
               "dimension of mu is incompatible")
  expect_error(dmulti_normal(1:3, mu = mu, Sigma = Sigma[1:2, 1:2],
                             check = TRUE),
               "dimension of Sigma is incompatible")
  expect_error(dmulti_normal(1:3, mu = mu, Sigma = Sigma[1:3, 3:1],
                             check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_normal(1.5, mu = mu, Sigma = Sigma, check = TRUE),
               "n must be a positive integer")
  expect_error(rmulti_normal(10, mu = mu, Sigma = Sigma[1:3, 3:1],
                             check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(dmulti_student(rnorm(3), mu = mu, Sigma = Sigma,
                              df = -1, check = TRUE),
               "df must be greater zero")
  expect_error(dmulti_student(rnorm(3), mu = mu, Sigma = Sigma[1:3, 3:1],
                              df = 30, check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_student(10, mu = mu, Sigma = Sigma,
                              df = -1, check = TRUE),
               "df must be greater zero")
})