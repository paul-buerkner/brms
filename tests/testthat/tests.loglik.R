test_that("loglik for location shift models works as expected", {
  ns <- 25
  s <- list(eta = matrix(rnorm(ns*2), ncol = 2),
            sigma = rchisq(ns, 3), nu = rgamma(ns, 4))
  data <- list(Y = rnorm(ns))
  
  ll_gaussian <- dnorm(x = data$Y[1], mean = s$eta[, 1], 
                       sd = s$sigma, log = TRUE)
  ll <- loglik_gaussian(1, data = data, samples = s)
  expect_equal(ll, ll_gaussian)
  
  ll_student <- dstudent(x = data$Y[2], df = s$nu, mu = 1 / s$eta[, 2], 
                         sigma = s$sigma, log = TRUE)
  ll <- loglik_student(2, data = data, samples = s, link = "inverse")
  expect_equal(ll, ll_student)
  
  ll_cauchy <- dstudent(x = data$Y[2], df = 1, mu = s$eta[, 2], 
                       sigma = s$sigma, log = TRUE)
  ll <- loglik_cauchy(2, data = data, samples = s)
  expect_equal(ll, ll_cauchy)
  
  # also test weighting
  data$weights <- sample(1:10, ns, replace = TRUE)
  ll <- loglik_gaussian(1, data = data, samples = s)
  expect_equal(ll, ll_gaussian * data$weights[1])
})

test_that("loglik for lognormal models works as expected", {
  ns <- 50
  s <- list(sigma = rchisq(ns, 3), eta = matrix(rnorm(ns*2), ncol = 2))
  data <- list(Y = rlnorm(ns))
  ll_lognormal <- dlnorm(x = data$Y[1], mean = s$eta[, 1], 
                         sd = s$sigma, log = TRUE)
  ll <- loglik_lognormal(1, data = data, samples = s)
  expect_equal(ll, ll_lognormal)
})

test_that("loglik for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), 
                dim = c(3, 3, 10))
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            Sigma = aperm(Sigma, c(3, 1, 2)), 
            nu = matrix(rgamma(ns, 5)))
  data <- list(Y = matrix(rnorm(nobs), ncol = nvars), 
               N_trait = ncols, K_trait = nvars)
  
  ll <- loglik_multi_gaussian(1, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_multi_student(2, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_multi_cauchy(2, data = data, samples = s)
  expect_equal(length(ll), ns)
})

test_that("loglik for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            sigma = matrix(rchisq(ns, 3)),
            nu = matrix(rgamma(ns, 5)),
            ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
            ma = matrix(rnorm(ns, 0.2, 1), ncol = 1))
  data <- list(Y = rnorm(nobs), begin_tg = 2, nrows_tg = 4,
               squared_se = rgamma(ns, 10))
  
  ll <- loglik_gaussian_cov(1, data = data, samples = s, link = "inverse")
  expect_equal(length(ll), ns)
  
  ll <- loglik_student_cov(1, data = data, samples = s, link = "log")
  expect_equal(length(ll), ns)
  
  ll <- loglik_cauchy_cov(1, data = data, samples = s)
  expect_equal(length(ll), ns)
})

test_that("loglik for count and survival models works correctly", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            shape = rgamma(ns, 4))
  data <- list(Y = rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1)), 
               max_obs = trials)
  i <- sample(1:nobs, 1)
  
  ll_binom <- dbinom(x = data$Y[i], prob = ilogit(s$eta[, i]), 
                     size = data$max_obs[i], log = TRUE)
  ll <- loglik_binomial(i, data = data, samples = s)
  expect_equal(ll, ll_binom)
  
  ll_pois <- dpois(x = data$Y[i], lambda = exp(s$eta[, i]), log = TRUE)
  ll <- loglik_poisson(i, data = data, samples = s)
  expect_equal(ll, ll_pois)
  
  ll_nbinom <- dnbinom(x = data$Y[i], mu = exp(s$eta[, i]), 
                       size = s$shape, log = TRUE)
  ll <- loglik_negbinomial(i, data = data, samples = s)
  expect_equal(ll, ll_nbinom)
  
  ll_geo <- dnbinom(x = data$Y[i], mu = exp(s$eta[, i]), 
                    size = 1, log = TRUE)
  ll <- loglik_geometric(i, data = data, samples = s)
  expect_equal(ll, ll_geo)
  
  ll_exp <- dexp(x = data$Y[i], rate = 1 / exp(s$eta[, i]), log = TRUE)
  ll <- loglik_exponential(i, data = data, samples = s)
  expect_equal(ll, ll_exp)
  
  ll_gamma <- dgamma(x = data$Y[i], shape = s$shape,
                     scale = exp(s$eta[, i]) / s$shape, log = TRUE)
  ll <- loglik_gamma(i, data = data, samples = s, link = "log")
  expect_equal(ll, ll_gamma)
  
  ll_weibull <- dweibull(x = data$Y[i], shape = s$shape,
                         scale = exp(s$eta[, i] / s$shape), log = TRUE)
  ll <- loglik_weibull(i, data = data, samples = s)
  expect_equal(ll, ll_weibull)
  
  ll_invgauss <- dinvgauss(x = data$Y[i], shape = s$shape,
                           mean = exp(s$eta[, i]), log = TRUE)
  ll <- loglik_inverse.gaussian(i, data = data, samples = s, link = "log")
  expect_equal(ll, ll_invgauss)
})

test_that("loglik for bernoulli and beta models works correctly", {
  ns <- 15
  nobs <- 10
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            phi = rgamma(ns, 4))
  i <- sample(1:nobs, 1)
  data <- list(Y = sample(0:1, nobs, replace = TRUE))
  ll_bern <- dbinom(x = data$Y[i], prob = ilogit(s$eta[, i]),
                    size = 1, log = TRUE)
  ll <- loglik_bernoulli(i, data = data, samples = s)
  expect_equal(ll, ll_bern)
  
  data <- list(Y = rbeta(nobs, 1, 1))
  ll_beta <- dbeta(x = data$Y[i], shape1 = ilogit(s$eta[, i]) * s$phi, 
                   shape2 = (1 - ilogit(s$eta[, i])) * s$phi, log = TRUE)
  ll <- loglik_beta(i, data = data, samples = s)
  expect_equal(ll, ll_beta)
})

test_that("loglik for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  resp <- rbinom(nobs-4, size = trials[1:4], prob = rbeta(nobs-4, 1, 1))
  s <- list(eta = matrix(rnorm(ns*nobs*2), ncol = nobs*2),
            shape = rgamma(ns, 4))
  data <- list(Y = c(resp, rep(0, 4)), N_trait = nobs, max_obs = trials)
  
  ll <- loglik_hurdle_poisson(1, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_hurdle_negbinomial(2, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_hurdle_gamma(5, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_zero_inflated_poisson(3, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_zero_inflated_binomial(4, data = data, samples = s)
  expect_equal(length(ll), ns)
  
  ll <- loglik_zero_inflated_negbinomial(6, data = data, samples = s)
  expect_equal(length(ll), ns)
})

test_that("loglik for categorical and ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 4
  s <- list(eta = array(rnorm(ns*nobs), dim = c(ns, nobs, ncat)))
  data <- list(Y = rep(1:ncat, 2), max_obs = ncat)
  ll <- sapply(1:nobs, loglik_categorical, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_cumulative, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_sratio, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_cratio, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_acat, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_acat, data = data, samples = s,
               link = "probit")
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("censored and truncated loglik run without errors", {
  ns <- 30
  nobs <- 3
  s <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
            sigma = rchisq(ns, 3))
  data <- list(Y = rnorm(ns), cens = c(-1,0,1))
  ll <- sapply(1:nobs, loglik_gaussian, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
  data <- list(Y = sample(-3:3, nobs), lb = -4, ub = 5)
  ll <- sapply(1:nobs, loglik_gaussian, data = data, samples = s)
  expect_equal(dim(ll), c(ns, nobs))
})