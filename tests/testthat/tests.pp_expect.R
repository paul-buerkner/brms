context("Tests for pp_expect helper functions")

# to reduce testing time on CRAN
skip_on_cran()

test_that("pp_expect helper functions run without errors", {
  # actually run pp_expect.brmsfit that call the helper functions
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  add_samples <- brms:::add_samples
  fit <- add_samples(fit, "shape", dist = "exp")
  fit <- add_samples(fit, "alpha", dist = "norm")
  fit <- add_samples(fit, "hu", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "zi", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "quantile", dist = "beta", shape1 = 2, shape2 = 1)
  fit <- add_samples(fit, "xi", dist = "unif", min = -1, max = 0.5)
  fit <- add_samples(fit, "ndt", dist = "exp")
  draws <- brms:::extract_draws(fit)
  draws$dpars$mu <- brms:::get_dpar(draws, "mu")
  draws$dpars$sigma <- brms:::get_dpar(draws, "sigma")
  draws$dpars$nu <- brms:::get_dpar(draws, "nu")
  nsamples <- nsamples(fit)
  nobs <- nobs(fit)
  
  # test preparation of truncated models
  draws$data$lb <- 0
  draws$data$ub <- 200
  mu <- brms:::pp_expect_trunc(draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # pseudo log-normal model
  fit$family <- fit$formula$family <- lognormal()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo shifted log-normal model
  fit$family <- fit$formula$family <- shifted_lognormal()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo skew-normal model
  fit$family <- fit$formula$family <- skew_normal()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo asym_laplace model
  fit$family <- fit$formula$family <- asym_laplace()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo zero_inflated_asym_laplace model
  fit$family <- fit$formula$family <- brmsfamily("zero_inflated_asym_laplace")
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo gen_extreme_value model
  fit$family <- fit$formula$family <- gen_extreme_value()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo weibull model
  fit$formula$pforms <- NULL
  fit$family <- fit$formula$family <- weibull()
  expect_equal(dim(SW(pp_expect(fit, summary = FALSE))), c(nsamples, nobs))
  
  # pseudo binomial model
  fit$autocor <- brms:::cor_empty()
  fit$family <- fit$formula$family <- binomial()
  expect_equal(dim(SW(pp_expect(fit, summary = FALSE))), c(nsamples, nobs))
  
  # pseudo hurdle poisson model
  fit$family <- fit$formula$family <- hurdle_poisson()
  fit$formula <- bf(count ~ Trt*Age + mo(Exp) + offset(Age) + (1+Trt|visit),
                    family = family(fit))
  expect_equal(dim(pp_expect(fit, summary = FALSE)), c(nsamples, nobs))
  
  # pseudo zero-inflated poisson model
  fit$family <- fit$formula$family <- zero_inflated_poisson()
  expect_equal(dim(pp_expect(fit, summary = FALSE)), c(nsamples, nobs))
  
  # pseudo custom model
  pp_expect_test <- function(draws) {
    draws$dpars$mu
  }
  fit$family <- fit$formula$family <- custom_family(
    "test", dpars = "mu", links = c("logit"),
    type = "int", vars = "trials[n]"
  )
  expect_equal(dim(pp_expect(fit, summary = FALSE)), c(nsamples, nobs))
  
  # truncated continuous models
  draws$dpars$shape <- c(as.matrix(fit, pars = "^shape$"))
  mu <- brms:::pp_expect_trunc_gaussian(draws, lb = 0, ub = 10)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::pp_expect_trunc_student(draws, lb = -Inf, ub = 15)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::pp_expect_trunc_lognormal(draws, lb = 2, ub = 15)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  draws$dpars$mu <- exp(draws$dpars$mu)
  mu <- brms:::pp_expect_trunc_gamma(draws, lb = 1, ub = 7)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::pp_expect_trunc_exponential(draws, lb = 0, ub = Inf)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- SW(brms:::pp_expect_trunc_weibull(draws, lb = -Inf, ub = Inf))
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # truncated discrete models
  data <- list(Y = sample(100, 10), trials = 1:10, N = 10)
  lb <- matrix(0, nrow = nsamples, ncol = nobs)
  ub <- matrix(100, nrow = nsamples, ncol = nobs)
  mu <- brms:::pp_expect_trunc_poisson(draws, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::pp_expect_trunc_negbinomial(draws, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::pp_expect_trunc_geometric(draws, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  draws$data$trials <- 120
  lb <- matrix(-Inf, nrow = nsamples, ncol = nobs)
  draws$dpars$mu <- brms:::ilink(draws$dpars$mu, "logit")
  mu <- brms:::pp_expect_trunc_binomial(draws, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
})

test_that("pp_expect_lagsar runs without errors", {
  draws <- list(
    dpars = list(mu = matrix(rnorm(30), nrow = 3)),
    ac = list(
      lagsar = matrix(c(0.3, 0.5, 0.7)), 
      Msar = matrix(1:100, 10, 10)
    ),
    nsamples = 3,
    nobs = 10,
    family = gaussian()
  )
  mu_new <- brms:::pp_expect_lagsar(draws)
  expect_equal(dim(mu_new), dim(draws$dpars$mu))
  expect_true(!identical(mu_new, draws$dpars$mu))
})

test_that("pp_expect for advanced count data distributions runs without errors", {
  ns <- 15
  nobs <- 5
  ncat <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu = array(rbeta(ns*nobs, 2, 2), dim = c(ns, nobs)),
    shape = array(rexp(ns*nobs, 3), dim = c(ns, nobs))
  )
  draws$family <- brmsfamily("discrete_weibull")
  pred <- suppressWarnings(brms:::pp_expect_discrete_weibull(draws))
  expect_equal(dim(pred), c(ns, nobs))
  
  draws$family <- brmsfamily("com_poisson")
  pred <- suppressWarnings(brms:::pp_expect_com_poisson(draws))
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("pp_expect for multinomial and dirichlet models runs without errors", {
  ns <- 15
  nobs <- 8
  ncat <- 3
  draws <- structure(list(nsamples = ns, nobs = nobs), class = "brmsdraws")
  draws$dpars <- list(
    mu1 = array(rnorm(ns*nobs), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs), dim = c(ns, nobs))
  )
  draws$data <- list(ncat = ncat, trials = sample(1:20, nobs))
 
  draws$family <- multinomial()
  pred <- brms:::pp_expect_multinomial(draws = draws)
  expect_equal(dim(pred), c(ns, nobs, ncat))
  
  draws$family <- dirichlet()
  pred <- brms:::pp_expect_dirichlet(draws = draws)
  expect_equal(dim(pred), c(ns, nobs, ncat))
})
