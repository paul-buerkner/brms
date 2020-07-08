context("Tests for distribution functions")

test_that("student distribution works correctly", {
  expect_equal(integrate(dstudent_t, -100, 100, df = 15, mu = 10, sigma = 5)$value, 1)
  expect_equal(dstudent_t(1, df = 10, mu = 0, sigma = 5), dt(1/5, df = 10)/5)
  expect_equal(pstudent_t(2, df = 20, mu = 2, sigma = 0.4), pt(0, df = 20))
  expect_equal(qstudent_t(0.7, df = 5, mu = 2, sigma = 3), 2 + 3*qt(0.7, df = 5))
  expect_equal(length(rstudent_t(10, df = 10, mu = rnorm(10), sigma = 1:10)), 10)
})

test_that("multivariate normal and student distributions work correctly", {
  mu <- rnorm(3)
  Sigma <- cov(matrix(rnorm(300), ncol = 3))
  expect_equal(dmulti_normal(1:3, mu = mu, Sigma = Sigma),
               mnormt::dmnorm(1:3, mu, Sigma))
  expect_equal(dmulti_student_t(1:3, mu = mu, Sigma = Sigma, df = 10, log = TRUE),
               mnormt::dmt(1:3, df = 10, mean = mu, S = Sigma, log = TRUE))
  expect_equal(dim(rmulti_normal(7, mu = mu, Sigma = Sigma)), c(7, 3))
  expect_equal(dim(rmulti_student_t(7, mu = mu, Sigma = Sigma, df = 10)), 
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
  expect_error(dmulti_student_t(rnorm(3), mu = mu, Sigma = Sigma,
                                df = -1, check = TRUE),
               "df must be greater than 0")
  expect_error(dmulti_student_t(rnorm(3), mu = mu, Sigma = Sigma[1:3, 3:1],
                                df = 30, check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_student_t(10, mu = mu, Sigma = Sigma,
                                df = -1, check = TRUE),
               "df must be greater than 0")
})

test_that("von_mises distribution functions run without errors", {
  n <- 10
  res <- dvon_mises(runif(n, -pi, pi), mu = 1, kappa = 1:n)
  expect_true(length(res) == n)
  res <- pvon_mises(runif(n, -pi, pi), mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
  res <- rvon_mises(n, mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
})

test_that("skew_normal distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dskew_normal(x, mu = 1, sigma = 2, alpha = 1)
  expect_true(length(res) == n)
  res <- pskew_normal(x, mu = rnorm(n), sigma = 1:n, 
                             alpha = 3, log.p = TRUE)
  expect_true(length(res) == n)
  res <- qskew_normal(x, mu = rnorm(n), sigma = 1:n, 
                      alpha = 3, log.p = TRUE)
  expect_true(length(res) == n)
  res <- rskew_normal(n, mu = rnorm(n), sigma = 10, alpha = -4:5)
  expect_true(length(res) == n)
})

test_that("exgaussian distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dexgaussian(x, mu = 1, sigma = 2, beta = 1)
  expect_true(length(res) == n)
  res <- pexgaussian(x, mu = rnorm(n), sigma = 1:n, 
                            beta = 3, log.p = TRUE)
  expect_true(length(res) == n)
  res <- rexgaussian(n, mu = rnorm(n), sigma = 10, beta = 1:10)
  expect_true(length(res) == n)
})

test_that("frechet distribution functions run without errors", {
  n <- 10
  x <- 21:30
  res <- dfrechet(x, loc = 1, scale = 2, shape = 1, log = TRUE)
  expect_true(length(res) == n)
  loc <- 1:10
  res <- pfrechet(x, loc = loc, scale = 1:n, shape = 3)
  expect_true(length(res) == n)
  q <- qfrechet(res, loc = loc, scale = 1:n, shape = 3)
  expect_equal(x, q)
  res <- rfrechet(n, loc = loc, scale = 10, shape = 1:10)
  expect_true(length(res) == n)
})

test_that("inv_gaussian distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- dinv_gaussian(x, mu = 1, shape = 1)
  expect_true(length(res) == n)
  res <- pinv_gaussian(x, mu = abs(rnorm(n)), shape = 3)
  expect_true(length(res) == n)
  res <- rinv_gaussian(n, mu = abs(rnorm(n)), shape = 1:10)
  expect_true(length(res) == n)
})

test_that("gen_extreme_value distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- dgen_extreme_value(x, mu = 1, sigma = 2, xi = 1)
  expect_true(length(res) == n)
  res <- pgen_extreme_value(x, mu = rnorm(n), sigma = 1:n, xi = 3)
  expect_true(length(res) == n)
  res <- rgen_extreme_value(n, mu = rnorm(n), sigma = 10, xi = 1:10)
  expect_true(length(res) == n)
})

test_that("asym_laplace distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dasym_laplace(x, mu = 1, sigma = 2, quantile = 0.5)
  expect_true(length(res) == n)
  res <- pasym_laplace(x, mu = rnorm(n), sigma = 1:n, quantile = 0.3)
  expect_true(length(res) == n)
  res <- rasym_laplace(n, mu = rnorm(n), sigma = 10, 
                              quantile = runif(n, 0, 1))
  expect_true(length(res) == n)
})

test_that("zero-inflated distribution functions run without errors", {
  n <- 10
  x <- rpois(n, lambda = 1)
  
  res <- dzero_inflated_poisson(x, lambda = 1, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_poisson(x, lambda = 1, zi = 0.1)
  expect_true(length(res) == n)
  
  res <- dzero_inflated_negbinomial(x, mu = 2, shape = 5, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_negbinomial(x, mu = 2, shape = 5, zi = 0.1)
  expect_true(length(res) == n)
  
  res <- dzero_inflated_binomial(x, size = c(2, 10), prob = 0.4, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_binomial(x, size = c(2, 10), prob = 0.4, zi = 0.1)
  expect_true(length(res) == n)
  
  x <- c(rbeta(n - 2, shape1 = 2, shape2 = 3), 0, 0)
  res <- dzero_inflated_beta(x, shape1 = 2, shape2 = 3, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_beta(x, shape1 = 2, shape2 = 3, zi = 0.1)
  expect_true(length(res) == n)
})

test_that("hurdle distribution functions run without errors", {
  n <- 10
  x <- rpois(n, lambda = 1)
  
  res <- dhurdle_poisson(x, lambda = 1, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_poisson(x, lambda = 1, hu = 0.1)
  expect_true(length(res) == n)
  
  res <- dhurdle_negbinomial(x, mu = 2, shape = 5, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_negbinomial(x, mu = 2, shape = 5, hu = 0.1)
  expect_true(length(res) == n)
  
  res <- dhurdle_gamma(x, shape = 1, scale = 3, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_gamma(x, shape = 1, scale = 3, hu = 0.1)
  expect_true(length(res) == n)
  
  res <- dhurdle_lognormal(x, mu = 2, sigma = 5, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_lognormal(x, mu = 2, sigma = 5, hu = 0.1)
  expect_true(length(res) == n)
})

test_that("wiener distribution functions run without errors", {
  set.seed(1234)
  n <- 10
  x <- seq(0.1, 1, 0.1)
  alpha <- rexp(n)
  tau <- 0.05
  beta <- 0.5
  delta <- rnorm(n)
  resp <- sample(c(0, 1), n, TRUE)
  
  d1 <- dwiener(x, alpha, tau, beta, delta, resp, backend = "Rwiener")
  d2 <- dwiener(x, alpha, tau, beta, delta, resp, backend = "rtdists")
  expect_equal(d1, d2)
  
  r1 <- rwiener(n, alpha, tau, beta, delta, backend = "Rwiener")
  r2 <- rwiener(n, alpha, tau, beta, delta, backend = "rtdists")
  expect_equal(names(r1), names(r2))
  expect_equal(dim(r1), dim(r2))
})
