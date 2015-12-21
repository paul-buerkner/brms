test_that("Test that loglik for (weighted) linear models works as expected", {
  ns <- 200
  s <- list(eta = matrix(rnorm(ns*2), ncol = 2),
            sigma = rchisq(ns, 3), nu = rgamma(ns, 4))
  data <- list(Y = rnorm(ns))
  ll_gaussian <- dnorm(x = data$Y[1], mean = s$eta[, 1], 
                       sd = s$sigma, log = TRUE)
  expect_equal(loglik_gaussian(1, data = data, samples = s), ll_gaussian)
  ll_student <- dstudent(x = data$Y[2], df = s$nu, mu = 1 / s$eta[, 2], 
                         sigma = s$sigma, log = TRUE)
  expect_equal(loglik_student(2, data = data, samples = s, link = "inverse"), 
               ll_student)
  ll_cauchy <- dstudent(x = data$Y[2], df = 1, mu = s$eta[, 2], 
                       sigma = s$sigma, log = TRUE)
  expect_equal(loglik_cauchy(2, data = data, samples = s), ll_cauchy)
  # also test weighting
  data$weights <- sample(1:10, ns, replace = TRUE)
  expect_equal(loglik_gaussian(1, data = data, samples = s), 
               ll_gaussian * data$weights[1])
})

test_that("Test that loglik for lognormal models works as expected", {
  ns <- 50
  s <- list(sigma = rchisq(ns, 3), eta = matrix(rnorm(ns*2), ncol = 2))
  data <- list(Y = rlnorm(ns))
  ll_lognormal <- dlnorm(x = data$Y[1], mean = s$eta[, 1], 
                         sd = s$sigma, log = TRUE)
  expect_equal(loglik_lognormal(1, data = data, samples = s), ll_lognormal)
})

test_that("Test that loglik for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(cbind(rnorm(100), rnorm(100), rnorm(100))), 
                dim = c(3, 3, 10))
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            Sigma = aperm(Sigma, c(3, 1, 2)), 
            nu = matrix(rgamma(ns, 5)))
  data <- list(Y = matrix(rnorm(nobs), ncol = nvars), 
               N_trait = ncols, K_trait = nvars)
  expect_equal(length(loglik_multi_gaussian(1, data = data, samples = s)), ns)
  expect_equal(length(loglik_multi_student(2, data = data, samples = s)), ns)
  expect_equal(length(loglik_multi_cauchy(2, data = data, samples = s)), ns)
})

test_that("Test that loglik for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            sigma = matrix(rchisq(ns, 3)),
            nu = matrix(rgamma(ns, 5)),
            ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
            ma = matrix(rnorm(ns, 0.2, 1), ncol = 1))
  data <- list(Y = rnorm(nobs), begin_tg = 2, nrows_tg = 4,
               squared_se = rgamma(ns, 10))
  expect_equal(length(loglik_gaussian_cov(1, data = data, samples = s,
                                          link = "inverse")), ns)
  expect_equal(length(loglik_student_cov(1, data = data, samples = s,
                                         link = "log")), ns)
  expect_equal(length(loglik_cauchy_cov(1, data = data, samples = s,
                                        link = "identity")), ns)
})

test_that("Test that loglik for count/survival models works correctly", {
  ns <- 200
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  s <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            shape = rgamma(ns, 4))
  data <- list(Y = rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1)), 
               max_obs = trials)
  i <- sample(1:nobs, 1)
  ll_binom <- dbinom(x = data$Y[i], prob = ilogit(s$eta[, i]), 
                     size = data$max_obs[i], log = TRUE)
  expect_equal(loglik_binomial(i, data = data, samples = s), ll_binom)
  ll_pois <- dpois(x = data$Y[i], lambda = exp(s$eta[, i]), log = TRUE)
  expect_equal(loglik_poisson(i, data = data, samples = s), ll_pois)
  ll_nbinom <- dnbinom(x = data$Y[i], mu = exp(s$eta[, i]), 
                       size = s$shape, log = TRUE)
  expect_equal(loglik_negbinomial(i, data = data, samples = s), ll_nbinom)
  ll_geo <- dnbinom(x = data$Y[i], mu = exp(s$eta[, i]), 
                    size = 1, log = TRUE)
  expect_equal(loglik_geometric(i, data = data, samples = s), ll_geo)
  ll_exp <- dexp(x = data$Y[i], rate = 1 / exp(s$eta[, i]), log = TRUE)
  expect_equal(loglik_exponential(i, data = data, samples = s), ll_exp)
  ll_gamma <- dgamma(x = data$Y[i], shape = s$shape,
                     scale = (1 / s$eta[, i]) / s$shape, log = TRUE)
  expect_equal(loglik_gamma(i, data = data, samples = s), ll_gamma)
  ll_weibull <- dweibull(x = data$Y[i], shape = s$shape,
                         scale = exp(s$eta[, i] / s$shape), log = TRUE)
  expect_equal(loglik_weibull(i, data = data, samples = s), ll_weibull)
  ll_invgauss <- dinvgauss(x = data$Y[i], shape = s$shape,
                           mean = exp(s$eta[, i]), log = TRUE)
  expect_equal(loglik_inverse.gaussian(i, data = data, samples = s,
                                       link = "log"), 
               ll_invgauss)
})

