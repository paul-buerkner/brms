context("Tests for posterior_predict helper functions")

test_that("posterior_predict for location shift models runs without errors", {
  ns <- 30
  nobs <- 10
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3), nu = rgamma(ns, 4)
  )
  i <- sample(nobs, 1)

  pred <- brms:::posterior_predict_gaussian(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_student(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for various skewed models runs without errors", {
  ns <- 50
  nobs <- 2
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    alpha = rnorm(ns), ndt = 1
  )
  pred <- brms:::posterior_predict_lognormal(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_shifted_lognormal(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_exgaussian(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_skew_normal(1, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for aysm_laplace models runs without errors", {
  ns <- 50
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3),
    quantile = rbeta(ns, 2, 1),
    mu = matrix(rnorm(ns*2), ncol = 2),
    zi = rbeta(ns, 10, 10)
  )
  pred <- brms:::posterior_predict_asym_laplace(1, prep = prep)
  expect_equal(length(pred), ns)
  pred <- brms:::posterior_predict_zero_inflated_asym_laplace(1, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
  prep <- structure(list(ndraws = ns), class = "mvbrmsprep")
  prep$mvpars <- list(
    Mu = array(rnorm(ns*nobs*nvars), dim = c(ns, nobs, nvars)),
    Sigma = aperm(Sigma, c(3, 1, 2))
  )
  prep$dpars <- list(nu = rgamma(ns, 5))
  prep$data <- list(N = nobs, N_trait = ncols)

  pred <- brms:::posterior_predict_gaussian_mv(1, prep = prep)
  expect_equal(dim(pred), c(ns, nvars))

  pred <- brms:::posterior_predict_student_mv(2, prep = prep)
  expect_equal(dim(pred), c(ns, nvars))
})

test_that("posterior_predict for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns*nobs), ncol = nobs),
    sigma = rchisq(ns, 3),
    nu = rgamma(ns, 5)
  )
  prep$ac <- list(
    ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
    ma = matrix(rnorm(ns, 0.2, 1), ncol = 1),
    begin_tg = c(1, 5, 12), end_tg = c(4, 11, 15)
  )
  prep$data <- list(se = rgamma(ns, 10))

  prep$family$fun <- "gaussian_time"
  pred <- brms:::posterior_predict_gaussian_time(1, prep = prep)
  expect_equal(length(pred), ns * 4)

  prep$family$fun <- "student_time"
  pred <- brms:::posterior_predict_student_time(2, prep = prep)
  expect_equal(length(pred), ns * 7)
})

test_that("loglik for SAR models runs without errors", {
  ns = 3
  prep <- structure(list(ndraws = ns, nobs = 10), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(30), nrow = ns),
    nu = rep(2, ns),
    sigma = rep(10, ns)
  )
  prep$ac <- list(lagsar = matrix(c(0.3, 0.5, 0.7)), Msar = diag(10))

  pred <- brms:::posterior_predict_gaussian_lagsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::posterior_predict_student_lagsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))

  prep$ac$errorsar <- prep$ac$lagsar
  prep$ac$lagsar <- NULL
  pred <- brms:::posterior_predict_gaussian_errorsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::posterior_predict_student_errorsar(1, prep = prep)
  expect_equal(dim(pred), c(3, 10))
})

test_that("posterior_predict for FCOR models runs without errors", {
  ns <- 3
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(nobs * ns), nrow = ns),
    sigma = rep(1, ns), nu = rep(2, ns)
  )
  prep$ac <- list(Mfcor = diag(nobs))
  pred <- brms:::posterior_predict_gaussian_fcor(1, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
  pred <- brms:::posterior_predict_student_fcor(1, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for count and survival models runs without errors", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns * nobs), ncol = nobs),
    shape = rgamma(ns, 4), xi = 0, phi = rgamma(ns, 1)
  )
  prep$dpars$nu <- prep$dpars$sigma <- prep$dpars$shape + 1
  prep$data <- list(trials = trials)
  i <- sample(nobs, 1)

  prep$dpars$mu <- brms:::inv_cloglog(prep$dpars$eta)
  pred <- brms:::posterior_predict_binomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_beta_binomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_discrete_weibull(i, prep = prep)
  expect_equal(length(pred), ns)

  prep$dpars$mu <- exp(prep$dpars$eta)
  pred <- brms:::posterior_predict_poisson(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_negbinomial(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_negbinomial2(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_geometric(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_com_poisson(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_exponential(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_gamma(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_frechet(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_inverse.gaussian(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_gen_extreme_value(i, prep = prep)
  expect_equal(length(pred), ns)

  prep$family$link <- "log"
  pred <- brms:::posterior_predict_weibull(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for bernoulli and beta models works correctly", {
  ns <- 17
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs)),
    phi = rgamma(ns, 4)
  )
  i <- sample(1:nobs, 1)

  pred <- brms:::posterior_predict_bernoulli(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_beta(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for xbeta models works correctly", {
  skip_if_not_installed("betareg")

  ns <- 17
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs)),
    phi = rgamma(ns, 4),
    kappa = rexp(ns, 5)
  )
  i <- sample(1:nobs, 1)

  pred <- brms:::posterior_predict_xbeta(i, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_xbeta(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = 2 * atan(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    kappa = rgamma(ns, 4)
  )
  i <- sample(seq_len(nobs), 1)
  pred <- brms:::posterior_predict_von_mises(i, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
    shape = rgamma(ns, 4), phi = rgamma(ns, 1),
    zi = rbeta(ns, 1, 1), coi = rbeta(ns, 5, 7)
  )
  prep$dpars$hu <- prep$dpars$zoi <- prep$dpars$zi
  prep$data <- list(trials = trials)

  prep$dpars$mu <- exp(prep$dpars$eta)
  pred <- brms:::posterior_predict_hurdle_poisson(1, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_hurdle_negbinomial(2, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_hurdle_gamma(5, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_poisson(3, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_negbinomial(6, prep = prep)
  expect_equal(length(pred), ns)

  prep$dpars$mu <- brms:::inv_logit(prep$dpars$eta)
  pred <- brms:::posterior_predict_zero_inflated_binomial(4, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_beta_binomial(6, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_inflated_beta(8, prep = prep)
  expect_equal(length(pred), ns)

  pred <- brms:::posterior_predict_zero_one_inflated_beta(7, prep = prep)
  expect_equal(length(pred), ns)
})

test_that("posterior_predict for ordinal models runs without errors", {
  ns <- 50
  nobs <- 8
  nthres <- 3
  ncat <- nthres + 1
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = array(rnorm(ns * nobs), dim = c(ns, nobs)),
    disc = rexp(ns),
    hu = rbeta(ns, 1, 1)
  )
  prep$thres$thres <- array(0, dim = c(ns, nthres))
  prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  prep$family$link <- "logit"

  prep$family$family <- "cumulative"
  pred <- sapply(1:nobs, brms:::posterior_predict_cumulative, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "sratio"
  pred <- sapply(1:nobs, brms:::posterior_predict_sratio, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "cratio"
  pred <- sapply(1:nobs, brms:::posterior_predict_cratio, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "acat"
  pred <- sapply(1:nobs, brms:::posterior_predict_acat, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$link <- "probit"
  pred <- sapply(1:nobs, brms:::posterior_predict_acat, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$family$family <- "hurdle_cumulative"
  pred <- sapply(1:nobs, brms:::posterior_predict_hurdle_cumulative, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for categorical and related models runs without erros", {
  set.seed(1234)
  ns <- 50
  nobs <- 8
  ncat <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu1 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs))
  )
  prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  prep$family <- categorical()
  prep$refcat <- 1
  pred <- sapply(1:nobs, brms:::posterior_predict_categorical, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$data$trials <- sample(1:20, nobs)
  prep$family <- multinomial()
  pred <- brms:::posterior_predict_multinomial(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))

  prep$data$trials <- sample(1:20, nobs)
  prep$dpars$phi <- rexp(ns, 1)
  prep$family <- dirichlet_multinomial()
  pred <- brms:::posterior_predict_dirichlet_multinomial(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))

  prep$dpars$phi <- rexp(ns, 1)
  prep$family <- dirichlet()
  pred <- brms:::posterior_predict_dirichlet(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))

  prep$family <- brmsfamily("dirichlet2")
  prep$dpars$mu1 <- rexp(ns, 10)
  prep$dpars$mu2 <- rexp(ns, 10)
  prep$dpars$mu3 <- rexp(ns, 10)
  pred <- brms:::posterior_predict_dirichlet2(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))

  prep$family <- brmsfamily("logistic_normal")
  prep$dpars <- list(
    mu2 = rnorm(ns),
    mu3 = rnorm(ns),
    sigma2 = rexp(ns, 10),
    sigma3 = rexp(ns, 10)
  )
  prep$lncor <- rbeta(ns, 2, 1)
  pred <- brms:::posterior_predict_logistic_normal(i = sample(1:nobs, 1), prep = prep)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))
})

test_that("truncated posterior_predict run without errors", {
  ns <- 30
  nobs <- 15
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )
  prep$refcat <- 1

  prep$data <- list(lb = sample(-(4:7), nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_gaussian, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$dpars$mu <- exp(prep$dpars$mu)
  prep$data <- list(ub = sample(70:80, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_poisson, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))

  prep$data <- list(lb = rep(0, nobs), ub = sample(70:75, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::posterior_predict_poisson, prep = prep)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_predict for the wiener diffusion model runs without errors", {
  skip("skip as long as RWiener fails on R-devel for 3.6.0")
  ns <- 5
  nobs <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    bs = rchisq(ns, 3), ndt = rep(0.5, ns),
    bias = rbeta(ns, 1, 1)
  )
  prep$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
  i <- sample(1:nobs, 1)
  expect_equal(nrow(brms:::posterior_predict_wiener(i, prep)), ns)
})

test_that("posterior_predict_custom runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rbeta(ns * nobs * 2, 1, 1), ncol = nobs * 2)
  )
  prep$data <- list(trials = rep(1, nobs))
  prep$family <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]"
  )
  posterior_predict_beta_binomial2 <- function(i, prep) {
    mu <- prep$dpars$mu[, i]
    rbinom(prep$ndraws, size = prep$data$trials[i], prob = mu)
  }
  expect_equal(length(brms:::posterior_predict_custom(sample(1:nobs, 1), prep)), ns)
})

test_that("posterior_predict_gaussian runs with various 'output' values without error", {
  fit <- rename_pars(brms:::brmsfit_example3)
  prep <- brms::prepare_predictions(fit)
  model_fit <- fit$fit@sim
  S <- model_fit$chains * (model_fit$iter - model_fit$warmup)
  i <- 1

  # random draws from Gaussian
  rpred <- brms:::posterior_predict_gaussian(i, prep = prep, output = "random")
  expect_equal(length(rpred), S)

  # compute PIT values (q = prep$data$Y[i])
  PITs <- brms:::posterior_predict_gaussian(i, prep = prep, output = "probability")
  expect_equal(length(PITs), S)
  expect_true(all(PITs >= 0 & PITs <= 1))
  
  # compute cdf based on custom 'q'
  qpred <- brms:::posterior_predict_gaussian(i, q = 15, prep = prep, output = "probability")
  expect_equal(length(qpred), S)
  expect_false(all(PITs == qpred))
  expect_true(all(qpred >= 0 & qpred <= 1))
})

test_that("truncated posterior_predict_gaussian runs with various 'output' values without error", {
  skip_if_not_installed("truncnorm")
  set.seed(1335)
  ns <- 30
  nobs <- 15
  i <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )
  prep$data <- list(
    Y = rnorm(nobs),
    lb = replicate(nobs, 0),
    ub = replicate(nobs, 10)
  )
 
  mu <- brms:::get_dpar(prep, "mu", i = i)
  sigma <- brms:::get_dpar(prep, "sigma", i = i)
  sigma <- brms:::add_sigma_se(sigma, prep, i = i)

  # compute cdf for truncated distribution
  obs_trunc_PITs <- brms:::posterior_predict_gaussian(i, prep = prep, output = "probability")
  expected_PITs <- truncnorm::ptruncnorm(q = prep$data$Y[i], a = prep$data$lb[i],
    b = prep$data$ub[i], mean = mu, sd = sigma)
  expect_equal(obs_trunc_PITs, expected_PITs)

  # take random draws from a truncated distribution
  rpred <- brms:::posterior_predict_gaussian(i, prep = prep, output = "random")
  expect_true(all(rpred >= prep$data$lb[i] & rpred <= prep$data$ub[i]))
})

test_that("posterior_predict_student runs with various 'output' values without error", {
  set.seed(1334)
  ns <- 30
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3), 
    nu = rgamma(ns, 4)
  )
  prep$data <- list(Y = rstudent_t(nobs, df = 3))
  i <- 8

  # random draws from non-truncated t
  rpred <- brms:::posterior_predict_student(i, prep = prep, output = "random")
  expect_equal(length(rpred), ns)

  # compute PIT values (q = prep$data$Y[i])
  PITs <- brms:::posterior_predict_student(i, prep = prep, output = "probability")
  expect_equal(length(PITs), ns)
  expect_true(all(PITs >= 0 & PITs <= 1))
  
  # compute cdf based on custom 'q'
  qpred <- brms:::posterior_predict_student(i, q = 15, prep = prep, output = "probability")
  expect_equal(length(qpred), ns)
  expect_false(all(PITs == qpred))
  expect_true(all(qpred >= 0 & qpred <= 1))

  prep$data$lb <- replicate(nobs, 0)  
  prep$data$ub <- replicate(nobs, 30)

  # random draws from truncated t
  rpred <- brms:::posterior_predict_student(i, prep = prep, output = "random")
  expect_true(all(rpred >= prep$data$lb[i] & rpred <= prep$data$ub[i]))

  # compute PIT values for truncated t (q = prep$data$Y[i])
  PITs_trunc <- brms:::posterior_predict_student(i, prep = prep, output = "probability")
  expect_equal(length(PITs_trunc), ns)
  expect_false(all(PITs == PITs_trunc))
  
  # compute cdf for truncated t based on custom 'q'
  qpred_trunc <- brms:::posterior_predict_student(i, q = 15, prep = prep, output = "probability")
  expect_equal(length(qpred_trunc), ns)
  expect_false(all(qpred == qpred_trunc))
})

test_that("posterior_predict_binomial works for different 'output' values without error", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns * nobs), ncol = nobs),
    shape = rgamma(ns, 4), xi = 0, phi = rgamma(ns, 1)
  )
  prep$dpars$nu <- prep$dpars$sigma <- prep$dpars$shape + 1
  i <- 3

  prep$dpars$mu <- brms:::inv_cloglog(prep$dpars$eta)
  
  prep$data <- list(
    trialsb = trials,
    Y = rbinom(nobs, size = trials, prob = prep$dpars$mu)
  )
  # random draws from binomial
  pred <- brms:::posterior_predict_binomial(i, prep = prep, output = "random")
  expect_equal(length(pred), ns)

  # compute PIT values (q = prep$data$trials[i])
  PITs <- brms:::posterior_predict_binomial(i, prep = prep, output = "probability")
  expect_equal(length(PITs), ns)
  expect_true(all(PITs >= 0 & PITs <= 1))

  # compute PIT values for custom 'q' (e.g., q = 5)
  qpred <- brms:::posterior_predict_binomial(i, q = 5, prep = prep, output = "probability")
  expect_equal(length(qpred), ns)
  expect_true(all(qpred >= 0 & qpred <= 1))
  expect_false(all(PITs == qpred))
})


test_that("posterior_predict_poisson works for different 'output' values without error", {
  set.seed(1386)
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = exp(matrix(rnorm(ns * nobs), ncol = nobs))
  )
  prep$data <- list(
    Y = rpois(nobs, lambda = prep$dpars$mu)
  )
  i <- 4

  pred <- brms:::posterior_predict_poisson(i, prep = prep, output = "random")
  expect_equal(length(pred), ns)

  PITs <- brms:::posterior_predict_poisson(i, prep = prep, output = "probability")
  expect_equal(length(PITs), ns)
  expect_true(all(PITs >= 0 & PITs <= 1))

  # truncation interval [1, 6]
  prep$data$lb <- replicate(nobs, 1) 
  prep$data$ub <- replicate(nobs, 6)

  rpred_trunc <- brms:::posterior_predict_poisson(i, prep = prep, output = "random", ntrys = 1000)
  # check whether invalid draws were returned
  # in case of invalid draws, the corresponding draw is a double and not an integer
  # this implementation is not ideal when posterior_predict is used by developers outside brms
  # would be better to return NA for invalid draws, or to throw an error if ntrys is exceeded or so
  rpred_trunc <- brms:::check_discrete_trunc_bounds(rpred_trunc, prep$data$lb[i], prep$data$ub[i])
  expect_equal(length(rpred_trunc), ns)
  expect_true(all(rpred_trunc >= prep$data$lb[i] & rpred_trunc <= prep$data$ub[i]))

  PITs_trunc <- brms:::posterior_predict_poisson(i, prep = prep, output = "probability")
  expect_equal(length(PITs_trunc), ns)
  expect_true(all(PITs_trunc >= 0 & PITs_trunc <= 1))

  # analytical sanity check for the cdf of a zero-truncated Poisson distribution
  # truncation interval [1, Inf]
  .zero_trunc_pois <- function(x, lambda) {
    (ppois(x, lambda = lambda) - exp(-lambda)) / (1 - exp(-lambda))
  }
  prep$dpars <- list(mu = 5)
  prep$data$lb <- 1.
  prep$data$ub <- Inf
  prep$data$Y <- rpois(1, lambda = prep$dpars$mu)

  obs_cdf <- brms:::posterior_predict_poisson(1, prep = prep, output = "probability")
  expected_cdf <- .zero_trunc_pois(prep$data$Y[1], lambda = prep$dpars$mu)

  expect_equal(obs_cdf, expected_cdf)
})