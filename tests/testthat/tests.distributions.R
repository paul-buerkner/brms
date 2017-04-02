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
               dmvnorm(1:3, mean = mu, sigma = Sigma))
  expect_equal(dmulti_student(1:3, mu = mu, Sigma = Sigma, df = 10, log = TRUE),
               dmvt(1:3, df = 10, delta = mu, sigma = Sigma, log = TRUE))
  expect_equal(dim(rmulti_normal(7, mu = mu, Sigma = Sigma)), c(7, 3))
  expect_equal(dim(rmulti_student(7, mu = mu, Sigma = Sigma, df = 10)), 
               c(7, 3))
  # test errors
  expect_error(dmulti_normal(1:3, mu = rnorm(2), Sigma = Sigma, check = TRUE),
               "Dimension of mu is incorrect")
  expect_error(dmulti_normal(1:3, mu = mu, Sigma = Sigma[1:2, 1:2],
                             check = TRUE),
               "Dimension of Sigma is incorrect")
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
               "df must be greater than 0")
  expect_error(dmulti_student(rnorm(3), mu = mu, Sigma = Sigma[1:3, 3:1],
                              df = 30, check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_student(10, mu = mu, Sigma = Sigma,
                              df = -1, check = TRUE),
               "df must be greater than 0")
})

test_that("von_mises distribution functions run without errors", {
  n <- 10
  res <- brms:::dvon_mises(runif(n, -pi, pi), mu = 1, kappa = 1:n)
  expect_true(length(res) == n)
  res <- brms:::pvon_mises(runif(n, -pi, pi), mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
  res <- brms:::rvon_mises(n, mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
})

test_that("exgaussian distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- brms:::dexgaussian(x, mu = 1, sigma = 2, beta = 1)
  expect_true(length(res) == n)
  res <- brms:::pexgaussian(x, mu = rnorm(n), sigma = 1:n, 
                            beta = 3, log.p = TRUE)
  expect_true(length(res) == n)
  res <- brms:::rexgaussian(n, mu = rnorm(n), sigma = 10, beta = 1:10)
  expect_true(length(res) == n)
})

test_that("frechet distribution functions run without errors", {
  n <- 10
  x <- 21:30
  res <- brms:::dfrechet(x, loc = 1, scale = 2, shape = 1, log = TRUE)
  expect_true(length(res) == n)
  loc <- 1:10
  res <- brms:::pfrechet(x, loc = loc, scale = 1:n, shape = 3)
  expect_true(length(res) == n)
  q <- brms:::qfrechet(res, loc = loc, scale = 1:n, shape = 3)
  expect_equal(x, q)
  res <- brms:::rfrechet(n, loc = loc, scale = 10, shape = 1:10)
  expect_true(length(res) == n)
})

test_that("inv_gaussian distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- brms:::dinv_gaussian(x, mu = 1, shape = 1)
  expect_true(length(res) == n)
  res <- brms:::pinv_gaussian(x, mu = abs(rnorm(n)), shape = 3)
  expect_true(length(res) == n)
  res <- brms:::rinv_gaussian(n, mu = abs(rnorm(n)), shape = 1:10)
  expect_true(length(res) == n)
})

test_that("gen_extreme_value distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- brms:::dgen_extreme_value(x, mu = 1, sigma = 2, xi = 1)
  expect_true(length(res) == n)
  res <- brms:::pgen_extreme_value(x, mu = rnorm(n), sigma = 1:n, xi = 3)
  expect_true(length(res) == n)
  res <- brms:::rgen_extreme_value(n, mu = rnorm(n), sigma = 10, xi = 1:10)
  expect_true(length(res) == n)
})

test_that("asym_laplace distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- brms:::dasym_laplace(x, mu = 1, sigma = 2, quantile = 0.5)
  expect_true(length(res) == n)
  res <- brms:::pasym_laplace(x, mu = rnorm(n), sigma = 1:n, quantile = 0.3)
  expect_true(length(res) == n)
  res <- brms:::rasym_laplace(n, mu = rnorm(n), sigma = 10, 
                              quantile = runif(n, 0, 1))
  expect_true(length(res) == n)
})
