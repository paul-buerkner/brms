context("Tests for posterior_epred helper functions")

# to reduce testing time on CRAN
skip_on_cran()

test_that("posterior_epred helper functions run without errors", {
  # actually run posterior_epred.brmsfit that call the helper functions
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  add_samples <- brms:::add_samples
  fit <- add_samples(fit, "shape", dist = "exp")
  fit <- add_samples(fit, "alpha", dist = "norm")
  fit <- add_samples(fit, "hu", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "zi", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "quantile", dist = "beta", shape1 = 2, shape2 = 1)
  fit <- add_samples(fit, "xi", dist = "unif", min = -1, max = 0.5)
  fit <- add_samples(fit, "ndt", dist = "exp")
  fit$formula$formula <- update(fit$formula$formula, .~. - arma(visit, patient))
  prep <- brms:::prepare_predictions(fit)
  prep$dpars$mu <- brms:::get_dpar(prep, "mu")
  prep$dpars$sigma <- brms:::get_dpar(prep, "sigma")
  prep$dpars$nu <- brms:::get_dpar(prep, "nu")
  nsamples <- nsamples(fit)
  nobs <- nobs(fit)
  
  # test preparation of truncated models
  prep$data$lb <- 0
  prep$data$ub <- 200
  mu <- brms:::posterior_epred_trunc(prep)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # pseudo log-normal model
  fit$family <- fit$formula$family <- lognormal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo shifted log-normal model
  fit$family <- fit$formula$family <- shifted_lognormal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo skew-normal model
  fit$family <- fit$formula$family <- skew_normal()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo asym_laplace model
  fit$family <- fit$formula$family <- asym_laplace()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo zero_inflated_asym_laplace model
  fit$family <- fit$formula$family <- brmsfamily("zero_inflated_asym_laplace")
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo gen_extreme_value model
  fit$family <- fit$formula$family <- gen_extreme_value()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo weibull model
  fit$formula$pforms <- NULL
  fit$family <- fit$formula$family <- weibull()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(nsamples, nobs))
  
  # pseudo binomial model
  fit$autocor <- brms:::cor_empty()
  fit$family <- fit$formula$family <- binomial()
  expect_equal(dim(SW(posterior_epred(fit, summary = FALSE))), c(nsamples, nobs))
  
  # pseudo hurdle poisson model
  fit$family <- fit$formula$family <- hurdle_poisson()
  fit$formula <- bf(count ~ Trt*Age + mo(Exp) + offset(Age) + (1+Trt|visit),
                    family = family(fit))
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(nsamples, nobs))
  
  # pseudo zero-inflated poisson model
  fit$family <- fit$formula$family <- zero_inflated_poisson()
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(nsamples, nobs))
  
  # pseudo custom model
  posterior_epred_test <- function(prep) {
    prep$dpars$mu
  }
  fit$family <- fit$formula$family <- custom_family(
    "test", dpars = "mu", links = c("logit"),
    type = "int", vars = "trials[n]"
  )
  expect_equal(dim(posterior_epred(fit, summary = FALSE)), c(nsamples, nobs))
  
  # truncated continuous models
  prep$dpars$shape <- c(as.matrix(fit, pars = "^shape$"))
  mu <- brms:::posterior_epred_trunc_gaussian(prep, lb = 0, ub = 10)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::posterior_epred_trunc_student(prep, lb = -Inf, ub = 15)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::posterior_epred_trunc_lognormal(prep, lb = 2, ub = 15)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  prep$dpars$mu <- exp(prep$dpars$mu)
  mu <- brms:::posterior_epred_trunc_gamma(prep, lb = 1, ub = 7)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::posterior_epred_trunc_exponential(prep, lb = 0, ub = Inf)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- SW(brms:::posterior_epred_trunc_weibull(prep, lb = -Inf, ub = Inf))
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # truncated discrete models
  data <- list(Y = sample(100, 10), trials = 1:10, N = 10)
  lb <- matrix(0, nrow = nsamples, ncol = nobs)
  ub <- matrix(100, nrow = nsamples, ncol = nobs)
  mu <- brms:::posterior_epred_trunc_poisson(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::posterior_epred_trunc_negbinomial(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- brms:::posterior_epred_trunc_geometric(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  prep$data$trials <- 120
  lb <- matrix(-Inf, nrow = nsamples, ncol = nobs)
  prep$dpars$mu <- brms:::ilink(prep$dpars$mu, "logit")
  mu <- brms:::posterior_epred_trunc_binomial(prep, lb = lb, ub = ub)
  expect_equal(dim(mu), c(nsamples, nobs))
})

test_that("posterior_epred_lagsar runs without errors", {
  prep <- list(
    dpars = list(mu = matrix(rnorm(30), nrow = 3)),
    ac = list(
      lagsar = matrix(c(0.3, 0.5, 0.7)), 
      Msar = matrix(1:100, 10, 10)
    ),
    nsamples = 3,
    nobs = 10,
    family = gaussian()
  )
  mu_new <- brms:::posterior_epred_lagsar(prep)
  expect_equal(dim(mu_new), dim(prep$dpars$mu))
  expect_true(!identical(mu_new, prep$dpars$mu))
})

test_that("posterior_epred for advanced count data distributions runs without errors", {
  ns <- 15
  nobs <- 5
  ncat <- 3
  prep <- structure(list(nsamples = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = array(rbeta(ns*nobs, 2, 2), dim = c(ns, nobs)),
    shape = array(rexp(ns*nobs, 3), dim = c(ns, nobs))
  )
  prep$family <- brmsfamily("discrete_weibull")
  pred <- suppressWarnings(brms:::posterior_epred_discrete_weibull(prep))
  expect_equal(dim(pred), c(ns, nobs))
  
  prep$family <- brmsfamily("com_poisson")
  pred <- suppressWarnings(brms:::posterior_epred_com_poisson(prep))
  expect_equal(dim(pred), c(ns, nobs))
})

test_that("posterior_epred for multinomial and dirichlet models runs without errors", {
  ns <- 15
  nobs <- 8
  ncat <- 3
  prep <- structure(list(nsamples = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu1 = array(rnorm(ns*nobs), dim = c(ns, nobs)),
    mu2 = array(rnorm(ns*nobs), dim = c(ns, nobs))
  )
  prep$data <- list(ncat = ncat, trials = sample(1:20, nobs))
 
  prep$family <- multinomial()
  pred <- brms:::posterior_epred_multinomial(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))
  
  prep$family <- dirichlet()
  pred <- brms:::posterior_epred_dirichlet(prep = prep)
  expect_equal(dim(pred), c(ns, nobs, ncat))
})
