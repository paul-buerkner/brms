context("Tests for predict helper functions")

test_that("predict for location shift models runs without errors", {
  ns <- 30
  nobs <- 10
  draws <- structure(list(nsamples = ns), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3), nu = rgamma(ns, 4)
  )
  i <- sample(nobs, 1)
  
  pred <- brms:::predict_gaussian(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_student(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for various skewed models runs without errors", {
  ns <- 50
  nobs <- 2
  draws <- structure(list(nsamples = ns), class = "brmsdraws")
  draws$dpars <- list(
    sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    alpha = rnorm(ns), ndt = 1
  )
  pred <- brms:::predict_lognormal(1, draws = draws)
  expect_equal(length(pred), ns)
  pred <- brms:::predict_shifted_lognormal(1, draws = draws)
  expect_equal(length(pred), ns)
  pred <- brms:::predict_exgaussian(1, draws = draws)
  expect_equal(length(pred), ns)
  pred <- brms:::predict_skew_normal(1, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for aysm_laplace models runs without errors", {
  ns <- 50
  draws <- structure(list(nsamples = ns), class = "brmsdraws")
  draws$dpars <- list(
    sigma = rchisq(ns, 3), 
    quantile = rbeta(ns, 2, 1),
    mu = matrix(rnorm(ns*2), ncol = 2),
    zi = rbeta(ns, 10, 10)
  )
  pred <- brms:::predict_asym_laplace(1, draws = draws)
  expect_equal(length(pred), ns)
  pred <- brms:::predict_zero_inflated_asym_laplace(1, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
  draws <- structure(list(nsamples = ns), class = "mvbrmsdraws")
  draws$mvpars <- list(
    Mu = array(rnorm(ns*nobs*nvars), dim = c(ns, nobs, nvars)),
    Sigma = aperm(Sigma, c(3, 1, 2))
  )
  draws$dpars <- list(nu = rgamma(ns, 5))
  draws$data <- list(N = nobs, N_trait = ncols)
  
  pred <- brms:::predict_gaussian_mv(1, draws = draws)
  expect_equal(dim(pred), c(ns, nvars))
  
  pred <- brms:::predict_student_mv(2, draws = draws)
  expect_equal(dim(pred), c(ns, nvars))
})

test_that("predict for ARMA covariance models runs without errors", {
  ns <- 20
  nobs <- 15
  draws <- structure(list(nsamples = ns), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns*nobs), ncol = nobs),
    sigma = rchisq(ns, 3),
    nu = rgamma(ns, 5)
  )
  draws$ac <- list(             
    ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
    ma = matrix(rnorm(ns, 0.2, 1), ncol = 1),
    begin_tg = c(1, 5, 12), end_tg = c(4, 11, 15)
  )
  draws$data <- list(se = rgamma(ns, 10))
  
  draws$family$fun <- "gaussian_cov"
  pred <- brms:::predict_gaussian_cov(1, draws = draws)
  expect_equal(length(pred), ns * 4)
  
  draws$family$fun <- "student_cov"
  pred <- brms:::predict_student_cov(2, draws = draws)
  expect_equal(length(pred), ns * 7)
})

test_that("loglik for SAR models runs without errors", {
  ns = 3
  draws <- structure(list(nsamples = ns, nobs = 10), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(30), nrow = ns),
    nu = rep(2, ns),
    sigma = rep(10, ns)
  )
  draws$ac <- list(lagsar = matrix(c(0.3, 0.5, 0.7)), W = diag(10))
  
  pred <- brms:::predict_gaussian_lagsar(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::predict_student_lagsar(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  
  draws$ac$errorsar <- draws$ac$lagsar
  draws$ac$lagsar <- NULL
  pred <- brms:::predict_gaussian_errorsar(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::predict_student_errorsar(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
})

test_that("predict for 'cor_fixed' models runs without errors", {
  ns <- 3
  draws <- structure(list(nsamples = ns), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(30), nrow = ns),
    nu = rep(2, ns)
  )
  draws$ac <- list(V = diag(10))
  pred <- brms:::predict_gaussian_fixed(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
  pred <- brms:::predict_student_fixed(1, draws = draws)
  expect_equal(dim(pred), c(3, 10))
})

test_that("predict for count and survival models runs without errors", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    eta = matrix(rnorm(ns*nobs), ncol = nobs),
    shape = rgamma(ns, 4), xi = 0
  )
  draws$dpars$nu <- draws$dpars$sigma <- draws$dpars$shape + 1
  draws$data <- list(trials = trials)
  i <- sample(nobs, 1)
  
  draws$dpars$mu <- brms:::inv_cloglog(draws$dpars$eta)
  pred <- brms:::predict_binomial(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_discrete_weibull(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$dpars$mu <- exp(draws$dpars$eta)
  pred <- brms:::predict_poisson(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_negbinomial(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_geometric(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_com_poisson(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_exponential(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_gamma(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_frechet(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_inverse.gaussian(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_gen_extreme_value(i, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$family$link <- "log"
  pred <- brms:::predict_weibull(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for bernoulli and beta models works correctly", {
  ns <- 17
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = 2 * nobs)),
    phi = rgamma(ns, 4)
  )
  i <- sample(1:nobs, 1)
  
  pred <- brms:::predict_bernoulli(i, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_beta(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = 2 * atan(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    kappa = rgamma(ns, 4)
  )
  i <- sample(seq_len(nobs), 1)
  pred <- brms:::predict_von_mises(i, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    eta = matrix(rnorm(ns * nobs * 2), ncol = nobs * 2),
    shape = rgamma(ns, 4), phi = rgamma(ns, 1),
    zi = rbeta(ns, 1, 1), coi = rbeta(ns, 5, 7)
  )
  draws$dpars$hu <- draws$dpars$zoi <- draws$dpars$zi
  draws$data <- list(trials = trials)
  
  draws$dpars$mu <- exp(draws$dpars$eta)
  pred <- brms:::predict_hurdle_poisson(1, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_hurdle_negbinomial(2, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_hurdle_gamma(5, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_zero_inflated_poisson(3, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_zero_inflated_negbinomial(6, draws = draws)
  expect_equal(length(pred), ns)
  
  draws$dpars$mu <- brms:::inv_logit(draws$dpars$eta)
  pred <- brms:::predict_zero_inflated_binomial(4, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_zero_inflated_beta(8, draws = draws)
  expect_equal(length(pred), ns)
  
  pred <- brms:::predict_zero_one_inflated_beta(7, draws = draws)
  expect_equal(length(pred), ns)
})

test_that("predict for ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  nthres <- 3
  ncat <- nthres + 1
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = array(rnorm(ns * nobs), dim = c(ns, nobs)),
    disc = rexp(ns)
  )
  draws$thres$thres <- array(0, dim = c(ns, nthres))
  draws$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  draws$family$link <- "logit"
  
  draws$family$family <- "cumulative"
  pred <- sapply(1:nobs, brms:::predict_cumulative, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$family$family <- "sratio"
  pred <- sapply(1:nobs, brms:::predict_sratio, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$family$family <- "cratio"
  pred <- sapply(1:nobs, brms:::predict_cratio, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$family$family <- "acat"
  pred <- sapply(1:nobs, brms:::predict_acat, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$family$link <- "probit"
  pred <- sapply(1:nobs, brms:::predict_acat, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("predict for categorical and related models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu1 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs, 0, 0.1), dim = c(ns, nobs))
  )
  draws$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  draws$family <- categorical()
  pred <- sapply(1:nobs, brms:::predict_categorical, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$data$trials <- sample(1:20, nobs)
  draws$family <- multinomial()
  pred <- brms:::predict_multinomial(i = sample(1:nobs, 1), draws = draws)
  expect_equal(dim(pred), c(ns, ncat))
  
  draws$dpars$phi <- rexp(ns, 1)
  draws$family <- dirichlet()
  pred <- brms:::predict_dirichlet(i = sample(1:nobs, 1), draws = draws)
  expect_equal(dim(pred), c(ns, ncat))
  expect_equal(rowSums(pred), rep(1, nrow(pred)))
})

test_that("truncated predict run without errors", {
  ns <- 30
  nobs <- 15
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )

  draws$data <- list(lb = sample(-(4:7), nobs, TRUE))
  pred <- sapply(1:nobs, brms:::predict_gaussian, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$dpars$mu <- exp(draws$dpars$mu)
  draws$data <- list(ub = sample(70:80, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::predict_poisson, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$data <- list(lb = rep(0, nobs), ub = sample(70:75, nobs, TRUE))
  pred <- sapply(1:nobs, brms:::predict_poisson, draws = draws)
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("predict for the wiener diffusion model runs without errors", {
  skip("skip as long as RWiener fails on R-devel for 3.6.0")
  ns <- 5
  nobs <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    bs = rchisq(ns, 3), ndt = rep(0.5, ns),
    bias = rbeta(ns, 1, 1)
  )
  draws$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
  i <- sample(1:nobs, 1)
  expect_equal(nrow(brms:::predict_wiener(i, draws)), ns)
})

test_that("predict_custom runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rbeta(ns * nobs * 2, 1, 1), ncol = nobs * 2)
  )
  draws$data <- list(trials = rep(1, nobs))
  draws$family <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]"
  )
  predict_beta_binomial2 <- function(i, draws) {
    mu <- draws$dpars$mu[, i]
    rbinom(draws$nsamples, size = draws$data$trials[i], prob = mu)
  }
  expect_equal(length(brms:::predict_custom(sample(1:nobs, 1), draws)), ns)
})
