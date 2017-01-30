test_that("fitted helper functions run without errors", {
  # actually run fitted.brmsfit that call the helper functions
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  add_samples <- brms:::add_samples
  fit <- add_samples(fit, "shape", dist = "exp")
  fit <- add_samples(fit, "hu", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "zi", dist = "beta", shape1 = 1, shape2 = 1)
  fit <- add_samples(fit, "quantile", dist = "beta", shape1 = 2, shape2 = 1)
  fit <- add_samples(fit, "xi", dist = "unif", min = -1, max = 0.5)
  draws <- brms:::extract_draws(fit)
  eta <- brms:::linear_predictor(draws)
  nsamples <- nsamples(fit)
  nobs <- nobs(fit)
  
  # test preparation of truncated models
  draws$data$lb <- 0
  draws$data$ub <- 200
  mu <- brms:::fitted_trunc(eta, draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # pseudo log-normal model
  fit$family <- lognormal()
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo asym_laplace model
  fit$family <- asym_laplace()
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo gen_extreme_value model
  fit$family <- gen_extreme_value()
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo weibull model
  fit$formula$pforms <- NULL
  fit$family <- weibull()
  expect_equal(dim(SW(fitted(fit, summary = FALSE))), 
               c(nsamples, nobs))
  
  # pseudo binomial model
  fit$autocor <- cor_arma()
  fit$family <- attr(fit$formula, "family") <- binomial()
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo hurdle poisson model
  fit$family <- hurdle_poisson()
  fit$formula <- bf(count ~ Trt*Age + mono(Exp) + offset(Age) + (1+Trt|visit))
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # pseudo zero-inflated poisson model
  fit$family <- attr(fit$formula, "family") <- zero_inflated_poisson()
  expect_equal(dim(fitted(fit, summary = FALSE)), 
               c(nsamples, nobs))
  
  # truncated continuous models
  draws$nu <- c(posterior_samples(fit, pars = "^nu$", as.matrix = TRUE))
  draws$shape <- c(posterior_samples(fit, pars = "^shape$", as.matrix = TRUE))
  mu <- fitted_trunc_gaussian(eta, lb = 0, ub = 10, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- fitted_trunc_student(eta, lb = -Inf, ub = 15, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- fitted_trunc_lognormal(eta, lb = 2, ub = 15, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  exp_eta <- exp(eta)
  mu <- fitted_trunc_gamma(exp_eta, lb = 1, ub = 7, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- fitted_trunc_exponential(exp_eta, lb = 0, ub = Inf, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- SW(fitted_trunc_weibull(exp_eta, lb = -Inf, ub = Inf, draws = draws))
  expect_equal(dim(mu), c(nsamples, nobs))
  
  # truncated discrete models
  data <- list(Y = sample(100, 10), trials = 1:10, N = 10)
  lb <- matrix(0, nrow = nsamples, ncol = nobs)
  ub <- matrix(100, nrow = nsamples, ncol = nobs)
  mu <- fitted_trunc_poisson(exp_eta, lb = lb, ub = ub, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- fitted_trunc_negbinomial(exp_eta, lb = lb, ub = ub, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  mu <- fitted_trunc_geometric(exp_eta, lb = lb, ub = ub, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
  
  draws$data$trials <- 120
  lb <- matrix(-Inf, nrow = nsamples, ncol = nobs)
  expit_eta <- ilink(eta, "logit")
  mu <- fitted_trunc_binomial(expit_eta, lb = lb, ub = ub, draws = draws)
  expect_equal(dim(mu), c(nsamples, nobs))
})
