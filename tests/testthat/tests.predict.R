test_that("predict for location shift models runs without errors", {
  ns <- 30
  nobs <- 10
  s <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
            sigma = rchisq(ns, 3), nu = rgamma(ns, 4))
  data <- list()
  i <- sample(nobs, 1)
  
  pred <- predict_gaussian(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_student(i, data = data, samples = s, link = "log")
  expect_equal(length(pred), ns)
  
  pred <- predict_cauchy(i, data = data, samples = s, link = "inverse")
  expect_equal(length(pred), ns)
})

test_that("predict for lognormal models runs without errors", {
  ns <- 50
  nobs <- 2
  s <- list(sigma = rchisq(ns, 3), 
            eta = matrix(rnorm(ns * nobs), ncol = nobs))
  pred <- predict_lognormal(1, data = list(), samples = s)
  expect_equal(length(pred), ns)
})

test_that("predict for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), 
                dim = c(3, 3, 10))
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            Sigma = aperm(Sigma, c(3, 1, 2)), 
            nu = matrix(rgamma(ns, 5)))
  data <- list(N = nobs, N_trait = ncols)
  
  pred <- predict_multi_gaussian(1, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nvars))
  
  pred <- predict_multi_student(2, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nvars))
  
  pred <- predict_multi_cauchy(3, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nvars))
})

test_that("predict for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            sigma = matrix(rchisq(ns, 3)),
            nu = matrix(rgamma(ns, 5)),
            ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
            ma = matrix(rnorm(ns, 0.2, 1), ncol = 1))
  data <- list(begin_tg = c(1, 5, 12), nrows_tg = c(4, 7, 3),
               squared_se = rgamma(ns, 10))
  
  pred <- predict_gaussian_cov(1, data = data, samples = s, link = "inverse")
  expect_equal(length(pred), ns * 4)
  
  pred <- predict_student_cov(2, data = data, samples = s[-4], link = "log")
  expect_equal(length(pred), ns * 7)
  
  pred <- predict_cauchy_cov(3, data = data, samples = s[-5])
  expect_equal(length(pred), ns * 3)
})

test_that("predict for count and survival models runs without errors", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            shape = rgamma(ns, 4))
  data <- list(max_obs = trials)
  i <- sample(nobs, 1)
  
  pred <- predict_binomial(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_poisson(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_negbinomial(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_geometric(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_exponential(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_gamma(i, data = data, samples = s, link = "log")
  expect_equal(length(pred), ns)
  
  pred <- predict_weibull(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_inverse.gaussian(i, data = data, samples = s, 
                                   link = "log")
  expect_equal(length(pred), ns)
})

test_that("predict for bernoulli and beta models works correctly", {
  ns <- 17
  nobs <- 10
  s <- list(eta = matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs),
            phi = rgamma(ns, 4))
  i <- sample(1:nobs, 1)
  data <- list()
  
  pred <- predict_bernoulli(i, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_bernoulli(i, data = list(N_trait = nobs), samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_beta(i, data = data, samples = s)
  expect_equal(length(pred), ns)
})

test_that("predict for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  s <- list(eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
            shape = rgamma(ns, 4), phi = rgamma(ns, 1))
  data <- list(N_trait = nobs, max_obs = trials)
  
  pred <- predict_hurdle_poisson(1, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_hurdle_negbinomial(2, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_hurdle_gamma(5, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_poisson(3, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_binomial(4, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_negbinomial(6, data = data, samples = s)
  expect_equal(length(pred), ns)
  
  pred <- predict_zero_inflated_beta(8, data = data, samples = s)
  expect_equal(length(pred), ns)
})

test_that("predict for categorical and ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 4
  s <- list(eta = array(rnorm(ns*nobs), dim = c(ns, nobs, ncat)))
  data <- list(Y = rep(1:ncat, 2), max_obs = ncat)
  
  pred <- sapply(1:nobs, predict_categorical, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_cumulative, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_sratio, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_cratio, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_acat, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  pred <- sapply(1:nobs, predict_acat, data = data, samples = s,
               link = "probit")
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("truncated predict run without errors", {
  ns <- 30
  nobs <- 15
  s <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
            sigma = rchisq(ns, 3))
  
  data <- list(lb = -4)
  pred <- sapply(1:nobs, predict_gaussian, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  data <- list(ub = 70)
  pred <- sapply(1:nobs, predict_poisson, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
  
  data <- list(lb = 0, ub = 70)
  pred <- sapply(1:nobs, predict_poisson, data = data, samples = s)
  expect_equal(dim(pred), c(ns, nobs))
})