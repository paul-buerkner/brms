context("Tests for log_lik helper functions")

test_that("log_lik for location shift models works as expected", {
  ns <- 25
  prep <- structure(list(), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * 2), ncol = 2),
    sigma = rchisq(ns, 3), nu = rgamma(ns, 4)
  )
  prep$family <- gaussian()
  prep$family$fun <- "gaussian"
  prep$data <- list(Y = rnorm(ns))

  ll_gaussian <- dnorm(
    x = prep$data$Y[1], mean = prep$dpars$mu[, 1],
    sd = prep$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_gaussian(1, prep = prep)
  expect_equal(ll, ll_gaussian)

  ll_student <- dstudent_t(
    x = prep$data$Y[2], df = prep$dpars$nu,
    mu = prep$dpars$mu[, 2],
    sigma = prep$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_student(2, prep = prep)
  expect_equal(ll, ll_student)

  # also test weighting
  prep$data$weights <- sample(1:10, ns, replace = TRUE)
  ll <- brms:::log_lik_gaussian(1, prep = prep)
  expect_equal(ll, ll_gaussian * prep$data$weights[1])
})

test_that("log_lik for various skewed normal models works as expected", {
  ns <- 50
  prep <- structure(list(), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
    mu = matrix(rnorm(ns*2), ncol = 2),
    alpha = rnorm(ns), ndt = 1
  )
  prep$data <- list(Y = rlnorm(ns))

  ll_lognormal <- dlnorm(
    x = prep$data$Y[1], mean = prep$dpars$mu[, 1],
    sd = prep$dpars$sigma, log = TRUE
  )
  ll <- brms:::log_lik_lognormal(1, prep = prep)
  expect_equal(ll, ll_lognormal)

  ll_shifted_lognormal <- dshifted_lnorm(
    x = prep$data$Y[1], mean = prep$dpars$mu[, 1],
    sd = prep$dpars$sigma, shift = prep$dpars$ndt, log = TRUE
  )
  ll <- brms:::log_lik_shifted_lognormal(1, prep = prep)
  expect_equal(ll, ll_shifted_lognormal)

  ll_exgaussian <- dexgaussian(
    x = prep$data$Y[1], mu = prep$dpars$mu[, 1],
    sigma = prep$dpars$sigma, beta = prep$dpars$beta, log = TRUE
  )
  ll <- brms:::log_lik_exgaussian(1, prep = prep)
  expect_equal(ll, ll_exgaussian)

  ll_skew_normal <- dskew_normal(
    x = prep$data$Y[1], mu = prep$dpars$mu[, 1],
    sigma = prep$dpars$sigma, alpha = prep$dpars$alpha, log = TRUE
  )
  ll <- as.numeric(brms:::log_lik_skew_normal(1, prep = prep))
  expect_equal(ll, ll_skew_normal)
})

test_that("log_lik of aysm_laplace models runs without errors", {
  ns <- 50
  prep <- structure(list(), class = "brmsprep")
  prep$dpars <- list(
    sigma = rchisq(ns, 3),
    quantile = rbeta(ns, 2, 1),
    mu = matrix(rnorm(ns*2), ncol = 2),
    zi = rbeta(ns, 10, 10)
  )
  prep$data <- list(Y = brms:::rasym_laplace(ns))
  ll <- brms:::log_lik_asym_laplace(1, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_zero_inflated_asym_laplace(1, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik for multivariate linear models runs without errors", {
  ns <- 10
  nvars <- 3
  ncols <- 4
  nobs <- nvars * ncols
  prep <- structure(list(), class = "mvbrmsprep")
  Sigma = array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
  prep$mvpars <- list(
    Mu = array(rnorm(ns*nobs*nvars), dim = c(ns, nobs, nvars)),
    Sigma = aperm(Sigma, c(3, 1, 2))
  )
  prep$dpars <- list(nu = rgamma(ns, 5))
  prep$ndraws <- ns
  prep$data <- list(Y = matrix(rnorm(nobs), ncol = nvars))

  ll <- brms:::log_lik_gaussian_mv(1, prep = prep)
  expect_equal(length(ll), ns)
  ll <- brms:::log_lik_student_mv(2, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik for ARMA models runs without errors", {
  ns <- 20
  nobs <- 15
  prep <- structure(list(ndraws = ns), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns*nobs), ncol = nobs),
    sigma = rchisq(ns, 3),
    nu = rgamma(ns, 5) + 15
  )
  prep$ac <- list(
    ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
    ma = matrix(rbeta(ns, 0.2, 1), ncol = 1),
    begin_tg = 2, end_tg = 5
  )
  prep$data <- list(Y = rnorm(nobs), se = rgamma(ns, 10))

  prep$family$fun <- "gaussian_time"
  ll <- brms:::log_lik_gaussian_time(1, prep = prep)
  expect_equal(dim(ll), c(ns, 4))
  prep$family$fun <- "student_time"
  ll <- brms:::log_lik_student_time(1, prep = prep)
  expect_equal(dim(ll), c(ns, 4))
})

test_that("log_lik for SAR models runs without errors", {
  prep <- structure(list(ndraws = 3, nobs = 10), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(30), nrow = 3),
    nu = rep(10, 3),
    sigma = rep(10, 3)
  )
  prep$ac <-  list(
    lagsar = matrix(c(0.3, 0.5, 0.7)),
    Msar = diag(10)
  )
  prep$data <- list(Y = rnorm(10))

  ll <- brms:::log_lik_gaussian_lagsar(1, prep = prep)
  expect_equal(dim(ll), c(3, 10))
  ll <- brms:::log_lik_student_lagsar(1, prep = prep)
  expect_equal(dim(ll), c(3, 10))

  prep$ac$errorsar <- prep$ac$lagsar
  prep$ac$lagsar <- NULL
  ll <- brms:::log_lik_gaussian_errorsar(1, prep = prep)
  expect_equal(dim(ll), c(3, 10))
  ll <- brms:::log_lik_student_errorsar(1, prep = prep)
  expect_equal(dim(ll), c(3, 10))
})

test_that("log_lik for FCOR models runs without errors", {
  ns <- 3
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(nobs * ns), nrow = ns),
    sigma = rep(1, ns),
    nu = rep(10, ns)
  )
  prep$ac <- list(Mfcor = diag(nobs))
  prep$data$Y <- rnorm(nobs)
  ll <- brms:::log_lik_gaussian_fcor(1, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))
  ll <- brms:::log_lik_student_fcor(1, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("log_lik for count and survival models works correctly", {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns*nobs), ncol = nobs),
    shape = rgamma(ns, 4),
    xi = runif(ns, -1, 0.5),
    phi = rgamma(ns, 1)
  )
  prep$dpars$sigma <- 1 / prep$dpars$shape
  prep$dpars$nu <- prep$dpars$shape + 1
  prep$data <- list(
    Y = rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1)),
    trials = trials
  )
  i <- sample(nobs, 1)

  prep$dpars$mu <- brms:::inv_logit(prep$dpars$eta)
  ll_binom <- dbinom(
    x = prep$data$Y[i], prob = prep$dpars$mu[, i],
    size = prep$data$trials[i], log = TRUE
  )
  ll <- brms:::log_lik_binomial(i, prep = prep)
  expect_equal(ll, ll_binom)

  ll_beta_binom <- dbeta_binomial(
    x = prep$data$Y[i], size = prep$data$trials[i],
    mu = prep$dpars$mu[, i], phi = prep$dpars$phi, log = TRUE
  )
  ll <- brms:::log_lik_beta_binomial(i, prep = prep)
  expect_equal(ll, ll_beta_binom)

  # don't test the actual values as they will be -Inf for this data
  ll <- brms:::log_lik_discrete_weibull(i, prep = prep)
  expect_equal(length(ll), ns)

  prep$dpars$mu <- exp(prep$dpars$eta)
  ll_pois <- dpois(
    x = prep$data$Y[i], lambda = prep$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_poisson(i, prep = prep)
  expect_equal(ll, ll_pois)

  ll_nbinom <- dnbinom(
    x = prep$data$Y[i], mu = prep$dpars$mu[, i],
    size = prep$dpars$shape, log = TRUE
  )
  ll <- brms:::log_lik_negbinomial(i, prep = prep)
  expect_equal(ll, ll_nbinom)

  ll <- brms:::log_lik_negbinomial2(i, prep = prep)
  expect_equal(ll, ll_nbinom)

  ll_geo <- dnbinom(
    x = prep$data$Y[i], mu = prep$dpars$mu[, i],
    size = 1, log = TRUE
  )
  ll <- brms:::log_lik_geometric(i, prep = prep)
  expect_equal(ll, ll_geo)

  ll_com_pois <- brms:::dcom_poisson(
    x = prep$data$Y[i], mu = prep$dpars$mu[, i],
    shape = prep$dpars$shape, log = TRUE
  )
  ll <- brms:::log_lik_com_poisson(i, prep = prep)
  expect_equal(ll, ll_com_pois)

  ll_exp <- dexp(
    x = prep$data$Y[i], rate = 1 / prep$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_exponential(i, prep = prep)
  expect_equal(ll, ll_exp)

  ll_gamma <- dgamma(
    x = prep$data$Y[i], shape = prep$dpars$shape,
    scale = prep$dpars$mu[, i] / prep$dpars$shape,
    log = TRUE
  )
  ll <- brms:::log_lik_gamma(i, prep = prep)
  expect_equal(ll, ll_gamma)

  scale <- prep$dpars$mu[, i] / gamma(1 - 1 / prep$dpars$nu)
  ll_frechet <- dfrechet(
    x = prep$data$Y[i], shape = prep$dpars$nu,
    scale = scale, log = TRUE
  )
  ll <- brms:::log_lik_frechet(i, prep = prep)
  expect_equal(ll, ll_frechet)

  ll_invgauss <- dinv_gaussian(
    x = prep$data$Y[i], shape = prep$dpars$shape,
    mu = prep$dpars$mu[, i], log = TRUE
  )
  ll <- brms:::log_lik_inverse.gaussian(i, prep = prep)
  expect_equal(ll, ll_invgauss)

  ll_weibull <- dweibull(
    x = prep$data$Y[i], shape = prep$dpars$shape,
    scale = prep$dpars$mu[, i] / gamma(1 + 1 / prep$dpars$shape),
    log = TRUE
  )
  ll <- brms:::log_lik_weibull(i, prep = prep)
  expect_equal(ll, c(ll_weibull))

  # keep test at the end
  prep$family$link <- "identity"
  prep$data$Y[i] <- 0
  ll_gen_extreme_value <- SW(dgen_extreme_value(
    x = prep$data$Y[i], mu = prep$dpars$mu[, i],
    sigma = prep$dpars$sigma, xi = prep$dpars$xi, log = TRUE
  ))
  ll <- SW(brms:::log_lik_gen_extreme_value(i, prep = prep))
  expect_equal(ll, ll_gen_extreme_value)
})

test_that("log_lik for bernoulli and beta models works correctly", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = brms:::inv_logit(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    phi = rgamma(ns, 4)
  )
  prep$data <- list(Y = sample(0:1, nobs, replace = TRUE))

  i <- sample(1:nobs, 1)
  ll_bern <- dbinom(
    x = prep$data$Y[i], prob = prep$dpars$mu[, i], size = 1, log = TRUE
  )
  ll <- brms:::log_lik_bernoulli(i, prep = prep)
  expect_equal(ll, ll_bern)

  prep$data <- list(Y = rbeta(nobs, 1, 1))
  ll_beta <- dbeta(
    x = prep$data$Y[i], shape1 = prep$dpars$mu[, i] * prep$dpars$phi,
    shape2 = (1 - prep$dpars$mu[, i]) * prep$dpars$phi, log = TRUE
  )
  ll <- brms:::log_lik_beta(i, prep = prep)
  expect_equal(ll, ll_beta)
})

test_that("log_lik for circular models runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = 2 * atan(matrix(rnorm(ns * nobs * 2), ncol = nobs * 2)),
    kappa = rgamma(ns, 4)
  )
  prep$data <- list(Y = runif(nobs, -pi, pi))
  i <- sample(seq_len(nobs), 1)
  ll <- brms:::log_lik_von_mises(i, prep = prep)
  expect_equal(length(ll), ns)
  prep$data$cens <- sample(-1:1, nobs, TRUE)
  ll <- brms:::log_lik_von_mises(i, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik for zero-inflated and hurdle models runs without erros", {
  ns <- 50
  nobs <- 8
  trials <- sample(10:30, nobs, replace = TRUE)
  resp <- rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1))
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    eta = matrix(rnorm(ns*nobs), ncol = nobs),
    shape = rgamma(ns, 4),
    phi = rgamma(ns, 1),
    zi = rbeta(ns, 1, 1),
    coi = rbeta(ns, 5, 7)
  )
  prep$dpars$hu <- prep$dpars$zoi <- prep$dpars$zi
  prep$data <- list(Y = c(resp, rep(0, 4)), trials = trials)

  prep$dpars$mu <- exp(prep$dpars$eta)
  ll <- brms:::log_lik_hurdle_poisson(1, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_hurdle_negbinomial(5, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_hurdle_gamma(2, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_hurdle_gamma(8, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_zero_inflated_poisson(3, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_zero_inflated_negbinomial(6, prep = prep)
  expect_equal(length(ll), ns)

  prep$dpars$mu <- brms:::inv_logit(prep$dpars$eta)
  ll <- brms:::log_lik_zero_inflated_binomial(4, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_zero_inflated_beta_binomial(7, prep = prep)
  expect_equal(length(ll), ns)

  prep$data$Y[1:nobs] <- rbeta(nobs / 2, 0.5, 4)
  ll <- brms:::log_lik_zero_inflated_beta(6, prep = prep)
  expect_equal(length(ll), ns)

  ll <- brms:::log_lik_zero_one_inflated_beta(7, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik for ordinal models runs without erros", {
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

  ll <- sapply(1:nobs, brms:::log_lik_cumulative, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  ll <- sapply(1:nobs, brms:::log_lik_sratio, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  ll <- sapply(1:nobs, brms:::log_lik_cratio, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  ll <- sapply(1:nobs, brms:::log_lik_acat, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$family$link <- "probit"
  ll <- sapply(1:nobs, brms:::log_lik_acat, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  ll <- brms:::log_lik_hurdle_cumulative(3, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik for categorical and related models runs without erros", {
  ns <- 50
  nobs <- 8
  ncat <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu1 = array(rnorm(ns*nobs), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs), dim = c(ns, nobs))
  )
  prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
  prep$family <- categorical()
  prep$refcat <- 1
  ll <- sapply(1:nobs, brms:::log_lik_categorical, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$data$Y <- matrix(
    sample(1:20, nobs * ncat, TRUE),
    nrow = nobs, ncol = ncat
  )
  prep$data$trials <- sample(1:20, nobs)
  prep$family <- multinomial()
  ll <- sapply(1:nobs, brms:::log_lik_multinomial, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$data$trials <- sample(1:20, nobs)
  prep$dpars$phi <- rexp(ns, 10)
  prep$family <- dirichlet_multinomial()
  ll <- sapply(1:nobs, brms:::log_lik_dirichlet_multinomial, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$data$Y <- prep$data$Y / rowSums(prep$data$Y)
  prep$dpars$phi <- rexp(ns, 10)
  prep$family <- dirichlet()
  ll <- sapply(1:nobs, brms:::log_lik_dirichlet, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$family <- brmsfamily("dirichlet2")
  prep$dpars$mu1 <- rexp(ns, 10)
  prep$dpars$mu2 <- rexp(ns, 10)
  prep$dpars$mu3 <- rexp(ns, 10)
  ll <- sapply(1:nobs, brms:::log_lik_dirichlet2, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))

  prep$family <- brmsfamily("logistic_normal")
  prep$dpars <- list(
    mu2 = rnorm(ns),
    mu3 = rnorm(ns),
    sigma2 = rexp(ns, 10),
    sigma3 = rexp(ns, 10)
  )
  prep$lncor <- rbeta(ns, 2, 1)
  ll <- sapply(1:nobs, brms:::log_lik_logistic_normal, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("censored and truncated log_lik run without errors", {
  ns <- 30
  nobs <- 3
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rnorm(ns * nobs), ncol = nobs),
    sigma = rchisq(ns, 3)
  )
  prep$data <- list(Y = rnorm(ns), cens = c(-1,0,1))
  ll <- sapply(1:nobs, brms:::log_lik_gaussian, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))
  prep$data <- list(Y = sample(-3:3, nobs), lb = -4, ub = 5)
  ll <- sapply(1:nobs, brms:::log_lik_gaussian, prep = prep)
  expect_equal(dim(ll), c(ns, nobs))
})

test_that("log_lik for the wiener diffusion model runs without errors", {
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
  expect_equal(length(brms:::log_lik_wiener(i, prep)), ns)
})

test_that("log_lik for the xbeta model runs without errors", {
  ns <- 50
  nobs <- 8
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rbeta(ns * nobs, 1.2, 2.3), ncol = nobs),
    phi = rexp(ns, 0.01),
    kappa = rexp(ns, 2)
  )
  prep$data <- list(Y = rbeta(nobs, 2, 3))
  ll <- brms:::log_lik_xbeta(3, prep = prep)
  expect_equal(length(ll), ns)
})

test_that("log_lik_custom runs without errors", {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(rbeta(ns * nobs * 2, 1, 1), ncol = nobs * 2)
  )
  prep$data <- list(
    Y = sample(0:1, nobs, replace = TRUE),
    trials = rep(1, nobs)
  )
  prep$family <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]"
  )
  log_lik_beta_binomial2 <- function(i, prep) {
    mu <- prep$dpars$mu[, i]
    dbinom(prep$data$Y[i], size = prep$data$trials[i], prob = mu)
  }
  expect_equal(length(brms:::log_lik_custom(sample(1:nobs, 1), prep)), ns)
})
