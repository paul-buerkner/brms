context("Tests for log_lik helper functions")

test_that("log_lik for location shift models works as expected", {
  ns <- 25
  draws <- structure(list(), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns * 2), ncol = 2),
    sigma = rchisq(ns, 3), nu = rgamma(ns, 4)
  )
  draws$family <- gaussian()
  draws$family$fun <- "gaussian"
  draws$data <- list(Y = rnorm(ns))
  
  ll_gaussian <- dnorm(
    x = draws$data$Y[1], mean = draws$dpars$mu[, 1], 
    sd = draws$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_gaussian(1, draws = draws)
  expect_equal(ll, ll_gaussian)
  
  ll_student <- dstudent_t(
    x = draws$data$Y[2], df = draws$dpars$nu, 
    mu = draws$dpars$mu[, 2], 
    sigma = draws$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_student(2, draws = draws)
  expect_equal(ll, ll_student)
  
  # also test weighting
  draws$data$weights <- sample(1:10, ns, replace = TRUE)
  ll <- brms:::log_lik_gaussian(1, draws = draws)
  expect_equal(ll, ll_gaussian * draws$data$weights[1])
})

test_that("log_lik for various skewed normal models works as expected", {
  ns <- 50
  draws <- structure(list(), class = "brmsdraws")
  draws$dpars <- list(
    sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
    mu = matrix(rnorm(ns*2), ncol = 2),
    alpha = rnorm(ns), ndt = 1
  )
  draws$data <- list(Y = rlnorm(ns))
  
  ll_lognormal <- dlnorm(
    x = draws$data$Y[1], mean = draws$dpars$mu[, 1], 
    sd = draws$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_lognormal(1, draws = draws)
  expect_equal(ll, ll_lognormal)
  
  ll_shifted_lognormal <- dshifted_lnorm(
    x = draws$data$Y[1], mean = draws$dpars$mu[, 1], 
    sd = draws$dpars$sigma, shift = draws$dpars$ndt, log = TRUE
  )
  ll <- brms:::log_lik_shifted_lognormal(1, draws = draws)
  expect_equal(ll, ll_shifted_lognormal)
  
  ll_exgaussian <- dexgaussian(
    x = draws$data$Y[1], mu = draws$dpars$mu[, 1],
    sigma = draws$dpars$sigma, beta = draws$dpars$beta, log = TRUE
  )
  ll <- brms:::log_lik_exgaussian(1, draws = draws)
  expect_equal(ll, ll_exgaussian)
  
  ll_skew_normal <- dskew_normal(
    x = draws$data$Y[1], mu = draws$dpars$mu[, 1],
    sigma = draws$dpars$sigma, alpha = draws$dpars$alpha, log = TRUE
  )
  ll <- as.numeric(brms:::log_lik_skew_normal(1, draws = draws))
  expect_equal(ll, ll_skew_normal)
})

test_that("log_lik of aysm_laplace models runs without errors", {
  ns <- 50
  draws <- structure(list(), class = "brmsdraws")
  draws$dpars <- list(
    sigma = rchisq(ns, 3), 
    quantile = rbeta(ns, 2, 1),
    mu = matrix(rnorm(ns*2), ncol = 2)
  )
  draws$data <- list(Y = brms:::rasym_laplace(ns))
  ll <- brms:::log_lik_asym_laplace(1, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("log_lik for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  draws <- structure(list(), class = "mvbrmsdraws")
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
  draws$mvpars <- list(
    Mu = array(rnorm(ns*nobs*nvars), dim = c(ns, nobs, nvars)),
    Sigma = aperm(Sigma, c(3, 1, 2))
  ) 
  draws$dpars <- list(nu = rgamma(ns, 5))
  draws$nsamples <- ns
  draws$data <- list(Y = matrix(rnorm(nobs), ncol = nvars))
  
  ll <- brms:::log_lik_gaussian_mv(1, draws = draws)
  expect_equal(length(ll), ns)
  ll <- brms:::log_lik_student_mv(2, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("log_lik for ARMA covariance models runs without errors", {
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
    ma = matrix(rbeta(ns, 0.2, 1), ncol = 1),
    begin_tg = 2, nobs_tg = 4
  )
  draws$data <- list(Y = rnorm(nobs), se = rgamma(ns, 10))

  draws$family$fun <- "gaussian_cov"
  ll <- brms:::log_lik_gaussian_cov(1, draws = draws)
  expect_equal(dim(ll), c(ns, 4))
  # draws$family$fun <- "student_cov"
  # ll <- brms:::log_lik_student_cov(1, draws = draws)
  # expect_equal(length(ll), ns)
})

test_that("log_lik for SAR models runs without errors", {
  draws <- structure(list(nsamples = 3, nobs = 10), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(30), nrow = 3),
    nu = rep(2, 3),
    sigma = rep(10, 3)
  )
  draws$ac <-  list(
    lagsar = matrix(c(0.3, 0.5, 0.7)),
    W = diag(10)
  )
  draws$data <- list(Y = rnorm(10))

  ll <- brms:::log_lik_gaussian_lagsar(1, draws = draws)
  expect_equal(dim(ll), c(3, 10))
  # ll <- brms:::log_lik_student_lagsar(1, draws = draws)
  # expect_equal(length(ll), 3)

  draws$ac$errorsar <- draws$ac$lagsar
  draws$ac$lagsar <- NULL
  ll <- brms:::log_lik_gaussian_errorsar(1, draws = draws)
  expect_equal(dim(ll), c(3, 10))
  # ll <- brms:::log_lik_student_errorsar(1, draws = draws)
  # expect_equal(length(ll), 3)
})

test_that("log_lik for 'cor_fixed' models runs without errors", {
  draws <- structure(list(nsamples = 3), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(30), nrow = 3),
    nu = rep(2, 3)
  )
  draws$ac$V <- diag(10)
  draws$data$Y <- rnorm(10)
  ll <- brms:::log_lik_gaussian_fixed(1, draws = draws)
  expect_equal(dim(ll), c(3, 10))
  # ll <- brms:::log_lik_student_fixed(1, draws = draws)
  # expect_equal(length(ll), 3)
})

test_that("log_lik for count and survival models works correctly", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    eta = matrix(rnorm(ns*nobs), ncol = nobs),
    shape = rgamma(ns, 4),
    xi = runif(ns, -1, 0.5)
  )
  draws$dpars$nu <- draws$dpars$sigma <- draws$dpars$shape + 1
  draws$data <- list(
    Y = rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1)), 
    trials = trials
  )
  i <- sample(nobs, 1)
  
  draws$dpars$mu <- brms:::inv_logit(draws$dpars$eta)
  ll_binom <- dbinom(
    x = draws$data$Y[i], prob = draws$dpars$mu[, i], 
    size = draws$data$trials[i], log = TRUE
  )
  ll <- brms:::log_lik_binomial(i, draws = draws)
  expect_equal(ll, ll_binom)
  
  # don't test the actual values as they will be -Inf for this data
  ll <- brms:::log_lik_discrete_weibull(i, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$dpars$mu <- exp(draws$dpars$eta)
  ll_pois <- dpois(
    x = draws$data$Y[i], lambda = draws$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_poisson(i, draws = draws)
  expect_equal(ll, ll_pois)
  
  ll_nbinom <- dnbinom(
    x = draws$data$Y[i], mu = draws$dpars$mu[, i], 
    size = draws$dpars$shape, log = TRUE
  )
  ll <- brms:::log_lik_negbinomial(i, draws = draws)
  expect_equal(ll, ll_nbinom)
  
  ll_geo <- dnbinom(
    x = draws$data$Y[i], mu = draws$dpars$mu[, i], 
    size = 1, log = TRUE
  )
  ll <- brms:::log_lik_geometric(i, draws = draws)
  expect_equal(ll, ll_geo)
  
  ll_exp <- dexp(
    x = draws$data$Y[i], rate = 1 / draws$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_exponential(i, draws = draws)
  expect_equal(ll, ll_exp)
  
  ll_gamma <- dgamma(
    x = draws$data$Y[i], shape = draws$dpars$shape,
    scale = draws$dpars$mu[, i] / draws$dpars$shape, 
    log = TRUE
  )
  ll <- brms:::log_lik_gamma(i, draws = draws)
  expect_equal(ll, ll_gamma)
  
  scale <- draws$dpars$mu[, i] / gamma(1 - 1 / draws$dpars$nu)
  ll_frechet <- dfrechet(
    x = draws$data$Y[i], shape = draws$dpars$nu,
    scale = scale, log = TRUE
  )
  ll <- brms:::log_lik_frechet(i, draws = draws)
  expect_equal(ll, ll_frechet)
  
  ll_invgauss <- dinv_gaussian(
    x = draws$data$Y[i], shape = draws$dpars$shape,
    mu = draws$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_inverse.gaussian(i, draws = draws)
  expect_equal(ll, ll_invgauss)
  
  ll_weibull <- dweibull(
    x = draws$data$Y[i], shape = draws$dpars$shape,
    scale = draws$dpars$mu[, i] / gamma(1 + 1 / draws$dpars$shape),
    log = TRUE
  )
  ll <- brms:::log_lik_weibull(i, draws = draws)
  expect_equal(ll, c(ll_weibull))
  
  # keep test at the end
  draws$family$link <- "identity"
  draws$data$Y[i] <- 0
  ll_gen_extreme_value <- SW(dgen_extreme_value(
    x = draws$data$Y[i], mu = draws$dpars$mu[, i],
    sigma = draws$dpars$nu, xi = draws$dpars$xi, log = TRUE
  ))
  ll <- SW(brms:::log_lik_gen_extreme_value(i, draws = draws))
  expect_equal(ll, ll_gen_extreme_value)
})

test_that("log_lik for bernoulli and beta models works correctly", {
  ns <- 15
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    phi = rgamma(ns, 4)
  )
  draws$data <- list(Y = sample(0:1, nobs, replace = TRUE))
  
  i <- sample(1:nobs, 1)
  ll_bern <- dbinom(
    x = draws$data$Y[i], prob = draws$dpars$mu[, i], size = 1, log = TRUE
  )
  ll <- brms:::log_lik_bernoulli(i, draws = draws)
  expect_equal(ll, ll_bern)
  
  draws$data <- list(Y = rbeta(nobs, 1, 1))
  ll_beta <- dbeta(
    x = draws$data$Y[i], shape1 = draws$dpars$mu[, i] * draws$dpars$phi, 
    shape2 = (1 - draws$dpars$mu[, i]) * draws$dpars$phi, log = TRUE
  )
  ll <- brms:::log_lik_beta(i, draws = draws)
  expect_equal(ll, ll_beta)
})

test_that("log_lik for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = 2 * atan(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    kappa = rgamma(ns, 4)
  )
  draws$data <- list(Y = runif(nobs, -pi, pi))
  i <- sample(seq_len(nobs), 1)
  ll <- brms:::log_lik_von_mises(i, draws = draws)
  expect_equal(length(ll), ns)
  draws$data$cens <- sample(-1:1, nobs, TRUE)
  ll <- brms:::log_lik_von_mises(i, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("log_lik for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  resp <- rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1))
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    eta = matrix(rnorm(ns*nobs), ncol = nobs),
    shape = rgamma(ns, 4), 
    phi = rgamma(ns, 1),
    zi = rbeta(ns, 1, 1), 
  coi = rbeta(ns, 5, 7)
  )
  draws$dpars$hu <- draws$dpars$zoi <- draws$dpars$zi
  draws$data <- list(Y = c(resp, rep(0, 4)), trials = trials)
  
  draws$dpars$mu <- exp(draws$dpars$eta)
  ll <- brms:::log_lik_hurdle_poisson(1, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_hurdle_negbinomial(5, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_hurdle_gamma(2, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_hurdle_gamma(8, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_zero_inflated_poisson(3, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_zero_inflated_negbinomial(6, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$dpars$mu <- brms:::inv_logit(draws$dpars$eta)
  ll <- brms:::log_lik_zero_inflated_binomial(4, draws = draws)
  expect_equal(length(ll), ns)
  
  draws$data$Y[1:nobs] <- rbeta(nobs / 2, 0.5, 4)
  ll <- brms:::log_lik_zero_inflated_beta(6, draws = draws)
  expect_equal(length(ll), ns)
  
  ll <- brms:::log_lik_zero_one_inflated_beta(7, draws = draws)
  expect_equal(length(ll), ns)
})

test_that("log_lik for ordinal models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 4
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = array(rnorm(ns*nobs), dim = c(ns, nobs, ncat)),
    disc = rexp(ns)
  )
  draws$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  draws$family$link <- "logit"
  
  ll <- sapply(1:nobs, brms:::log_lik_cumulative, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, brms:::log_lik_sratio, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, brms:::log_lik_cratio, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  ll <- sapply(1:nobs, brms:::log_lik_acat, data = data, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  draws$family$link <- "probit"
  ll <- sapply(1:nobs, brms:::log_lik_acat, data = data, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("log_lik for categorical and related models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu1 = array(rnorm(ns*nobs), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs), dim = c(ns, nobs))
  )
  draws$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  draws$family <- categorical()
  ll <- sapply(1:nobs, brms:::log_lik_categorical, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  draws$data$Y <- matrix(
    sample(1:20, nobs * ncat, TRUE), 
    nrow = nobs, ncol = ncat
  )
  draws$data$trials <- sample(1:20, nobs)
  draws$family <- multinomial()
  ll <- sapply(1:nobs, brms:::log_lik_multinomial, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  
  draws$data$Y <- draws$data$Y / rowSums(draws$data$Y)
  draws$dpars$phi <- rexp(ns, 10)
  draws$family <- dirichlet()
  ll <- sapply(1:nobs, brms:::log_lik_dirichlet, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("censored and truncated log_lik run without errors", {
  ns <- 30
  nobs <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )
  draws$data <- list(Y = rnorm(ns), cens = c(-1,0,1))
  ll <- sapply(1:nobs, brms:::log_lik_gaussian, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
  draws$data <- list(Y = sample(-3:3, nobs), lb = -4, ub = 5)
  ll <- sapply(1:nobs, brms:::log_lik_gaussian, draws = draws)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("log_lik for the wiener diffusion model runs without errors", {
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
  expect_equal(length(brms:::log_lik_wiener(i, draws)), ns)
})

test_that("log_lik_custom runs without errors", {
  ns <- 15
  nobs <- 10
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = matrix(rbeta(ns * nobs * 2, 1, 1), ncol = nobs * 2)
  )
  draws$data <- list(
    Y = sample(0:1, nobs, replace = TRUE),
    trials = rep(1, nobs)
  )
  draws$family <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]"
  )
  log_lik_beta_binomial2 <- function(i, draws) {
    mu <- draws$dpars$mu[, i]
    dbinom(draws$data$Y[i], size = draws$data$trials[i], prob = mu)
  }
  expect_equal(length(brms:::log_lik_custom(sample(1:nobs, 1), draws)), ns)
})
