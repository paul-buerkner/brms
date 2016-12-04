test_that("loglik for location shift models works as expected", {
  ns <- 25
  draws <- list(eta = matrix(rnorm(ns * 2), ncol = 2),
                sigma = rchisq(ns, 3), nu = rgamma(ns, 4))
  draws$f$link <- "identity"
  draws$data <- list(Y = rnorm(ns))
  
  ll_gaussian <- dnorm(x = draws$data$Y[1], mean = draws$eta[, 1], 
                       sd = draws$sigma, log = TRUE)
  ll <- loglik_gaussian(1, draws = draws)
  expect_equal(ll, as.matrix(ll_gaussian))
  
  ll_cauchy <- dstudent(x = draws$data$Y[2], df = 1, mu = draws$eta[, 2], 
                        sigma = draws$sigma, log = TRUE)
  ll <- loglik_cauchy(2, draws = draws)
  expect_equal(ll, as.matrix(ll_cauchy))
  
  ll_student <- dstudent(x = draws$data$Y[2], df = draws$nu, 
                         mu = 1 / draws$eta[, 2], 
                         sigma = draws$sigma, log = TRUE)
  draws$f$link <- "inverse"
  ll <- loglik_student(2, draws = draws)
  expect_equal(ll, as.matrix(ll_student))
  
  # also test weighting
  draws$f$link <- "identity"
  draws$data$weights <- sample(1:10, ns, replace = TRUE)
  ll <- loglik_gaussian(1, draws = draws)
  expect_equal(ll, as.matrix(ll_gaussian * draws$data$weights[1]))
})

test_that("loglik for lognormal and exgaussian models works as expected", {
  ns <- 50
  draws <- list(sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
                eta = matrix(rnorm(ns*2), ncol = 2),
                f = lognormal())
  draws$data <- list(Y = rlnorm(ns))
  ll_lognormal <- dlnorm(x = draws$data$Y[1], mean = draws$eta[, 1], 
                         sd = draws$sigma, log = TRUE)
  ll <- loglik_lognormal(1, draws = draws)
  expect_equal(ll, as.matrix(ll_lognormal))
  
  ll_exgaussian <- dexgauss(x = draws$data$Y[1], mu = draws$eta[, 1], 
                            sigma = draws$sigma, beta = draws$beta,
                            log = TRUE)
  ll <- loglik_exgaussian(1, draws = draws)
  expect_equal(ll, ll_exgaussian)
})

test_that("loglik for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), 
                dim = c(3, 3, 10))
  draws <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
            Sigma = aperm(Sigma, c(3, 1, 2)), 
            nu = matrix(rgamma(ns, 5)),
            nsamples = ns)
  draws$data <- list(Y = matrix(rnorm(nobs), ncol = nvars), 
                     N_trait = ncols, K_trait = nvars)
  draws$f$link <- "identity"
  
  ll <- loglik_gaussian_mv(1, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_student_mv(2, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_cauchy_mv(2, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("loglik for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  draws <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
                sigma = matrix(rchisq(ns, 3)),
                nu = matrix(rgamma(ns, 5)),
                ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
                ma = matrix(rnorm(ns, 0.2, 1), ncol = 1),
                nsamples = ns)
  draws$data <- list(Y = rnorm(nobs), begin_tg = 2, nobs_tg = 4,
                     se2 = rgamma(ns, 10))
  
  draws$f$link <- "inverse"
  ll <- loglik_gaussian_cov(1, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$f$link <- "log"
  ll <- loglik_student_cov(1, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$f$link <- "identity"
  ll <- loglik_cauchy_cov(1, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("loglik for 'cor_fixed' models runs without errors", {
  draws <- list(eta = matrix(rnorm(30), nrow = 3),
                nu = matrix(rep(2, 3)),
                nsamples = 3)
  draws$data <- list(Y = rnorm(10), V = diag(10))
  draws$f$link <- "identity"
  ll <- loglik_gaussian_fixed(1, draws = draws)
  expect_equal(length(ll), 3)
  ll <- loglik_student_fixed(1, draws = draws)
  expect_equal(length(ll), 3)
  ll <- loglik_cauchy_fixed(1, draws = draws)
  expect_equal(length(ll), 3)
})

test_that("loglik for count and survival models works correctly", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- list(eta = matrix(rnorm(ns*nobs), ncol = nobs),
                shape = matrix(rgamma(ns, 4)), nsamples = ns)
  draws$data <- list(Y = rbinom(nobs, size = trials, 
                                prob = rbeta(nobs, 1, 1)), 
                     max_obs = trials)
  
  i <- sample(nobs, 1)
  
  draws$f$link <- "logit"
  ll_binom <- dbinom(x = draws$data$Y[i], prob = inv_logit(draws$eta[, i]), 
                     size = draws$data$max_obs[i], log = TRUE)
  ll <- loglik_binomial(i, draws = draws)
  expect_equal(ll, as.matrix(ll_binom))
  
  draws$f$link <- "log"
  ll_pois <- dpois(x = draws$data$Y[i], lambda = exp(draws$eta[, i]), 
                   log = TRUE)
  ll <- loglik_poisson(i, draws = draws)
  expect_equal(ll, as.matrix(ll_pois))
  
  ll_nbinom <- dnbinom(x = draws$data$Y[i], mu = exp(draws$eta[, i]), 
                       size = draws$shape, log = TRUE)
  ll <- loglik_negbinomial(i, draws = draws)
  expect_equal(ll, as.matrix(ll_nbinom))
  
  ll_geo <- dnbinom(x = draws$data$Y[i], mu = exp(draws$eta[, i]), 
                    size = 1, log = TRUE)
  ll <- loglik_geometric(i, draws = draws)
  expect_equal(ll, as.matrix(ll_geo))
  
  ll_exp <- dexp(x = draws$data$Y[i], rate = 1 / exp(draws$eta[, i]), 
                 log = TRUE)
  ll <- loglik_exponential(i, draws = draws)
  expect_equal(ll, as.matrix(ll_exp))
  
  ll_gamma <- dgamma(x = draws$data$Y[i], shape = draws$shape,
                     scale = exp(draws$eta[, i]) / draws$shape, log = TRUE)
  ll <- loglik_gamma(i, draws = draws)
  expect_equal(ll, as.matrix(ll_gamma))
  
  ll_weibull <- dweibull(x = draws$data$Y[i], shape = draws$shape,
                         scale = exp(draws$eta[, i] / draws$shape), log = TRUE)
  ll <- loglik_weibull(i, draws = draws)
  expect_equal(ll, as.matrix(ll_weibull))
  
  ll_invgauss <- dinvgauss(x = draws$data$Y[i], shape = draws$shape,
                           mean = exp(draws$eta[, i]), log = TRUE)
  ll <- loglik_inverse.gaussian(i, draws = draws)
  expect_equal(ll, ll_invgauss)
})

test_that("loglik for bernoulli and beta models works correctly", {
  ns <- 15
  nobs <- 10
  draws <- list(eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
                phi = matrix(rgamma(ns, 4)))
  draws$data <- list(Y = sample(0:1, nobs, replace = TRUE))
  draws$f$link <- "logit"
  i <- sample(1:nobs, 1)
  ll_bern <- dbinom(x = draws$data$Y[i], prob = inv_logit(draws$eta[, i]),
                    size = 1, log = TRUE)
  ll <- loglik_bernoulli(i, draws = draws)
  expect_equal(ll, as.matrix(ll_bern))
  
  draws$data <- list(Y = rbeta(nobs, 1, 1))
  ll_beta <- dbeta(x = draws$data$Y[i], shape1 = inv_logit(draws$eta[, i]) * draws$phi, 
                   shape2 = (1 - inv_logit(draws$eta[, i])) * draws$phi, log = TRUE)
  ll <- loglik_beta(i, draws = draws)
  expect_equal(ll, as.matrix(ll_beta))
})

test_that("loglik for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- list(eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
                kappa = matrix(rgamma(ns, 4)))
  draws$data <- list(Y = runif(nobs, -pi, pi))
  draws$f$link <- "tan_half"
  i <- sample(seq_len(nobs), 1)
  ll <- loglik_von_mises(i, draws = draws)
  expect_equal(length(ll), ns)
  draws$data$cens <- sample(-1:1, nobs, TRUE)
  ll <- loglik_von_mises(i, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("loglik for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  resp <- rbinom(nobs / 2, size = trials[1:(nobs / 2)], 
                 prob = rbeta(nobs / 2, 1, 1))
  draws <- list(eta = matrix(rnorm(ns*nobs*2), ncol = nobs*2),
                shape = matrix(rgamma(ns, 4)), 
                phi = matrix(rgamma(ns, 1)))
  draws$data <- list(Y = c(resp, rep(0, 4)), N_trait = nobs, 
                     max_obs = trials)
  draws$f$link <- "log"
  
  ll <- loglik_hurdle_poisson(1, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_hurdle_negbinomial(5, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_hurdle_gamma(2, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_hurdle_gamma(8, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_zero_inflated_poisson(3, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- loglik_zero_inflated_negbinomial(6, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$f$link <- "logit"
  ll <- loglik_zero_inflated_binomial(4, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$data$Y[1:(nobs / 2)] <- rbeta(nobs / 2, 0.5, 4)
  ll <- loglik_zero_inflated_beta(6, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("loglik for categorical and ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 4
  draws <- list(eta = array(rnorm(ns*nobs), dim = c(ns, nobs, ncat)),
                nsamples = ns)
  draws$data <- list(Y = rep(1:ncat, 2), max_obs = ncat)
  draws$f$link <- "logit"
  ll <- sapply(1:nobs, loglik_categorical, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_cumulative, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_sratio, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_cratio, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, loglik_acat, data = data, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  draws$f$link <- "probit"
  ll <- sapply(1:nobs, loglik_acat, data = data, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("censored and truncated loglik run without errors", {
  ns <- 30
  nobs <- 3
  draws <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
                sigma = matrix(rchisq(ns, 3)))
  draws$data <- list(Y = rnorm(ns), cens = c(-1,0,1))
  draws$f$link <- "identity"
  ll <- sapply(1:nobs, loglik_gaussian, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  draws$data <- list(Y = sample(-3:3, nobs), lb = -4, ub = 5)
  ll <- sapply(1:nobs, loglik_gaussian, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("loglik for the wiener diffusion model runs without errors", {
  ns <- 5
  nobs <- 3
  draws <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
                bs = matrix(rchisq(ns, 3)), ndt = matrix(rep(0.5, ns)),
                bias = matrix(rbeta(ns, 1, 1)))
  draws$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
  draws$f$link <- "identity"
  i <- sample(1:nobs, 1)
  expect_equal(length(loglik_wiener(i, draws)), ns)
})
