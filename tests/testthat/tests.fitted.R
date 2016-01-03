test_that("fitted helper functions run without errors", {
  # actually run fitted.brmsfit that call the helper functions
  fit <- rename_pars(brmsfit_example)
  fit <- add_samples(fit, "sigma_count", dist = "exp")
  fit <- add_samples(fit, "shape", dist = "exp")
  fit <- add_samples(fit, "nu", dist = "exp")
  eta <- linear_predictor(fit)
  # pseudo binomial model
  fit$family <- binomial()
  expect_equal(dim(fitted(fit, summary = FALSE)), c(80, 236))
  # pseudo log-normal model
  fit$family <- gaussian("log")
  expect_equal(dim(fitted(fit, summary = FALSE)), c(80, 236))
  # pseudo weibull model
  fit$family <- weibull()
  expect_equal(dim(fitted(fit, summary = FALSE)), c(80, 236))
  # pseudo hurdle poisson model
  names(fit$data)[1] <- "response"
  fit$family <- hurdle_poisson()
  expect_equal(dim(fitted(fit, summary = FALSE)), c(80, 118))
  # pseudo zero-inflated poisson model
  fit$family <- zero_inflated_poisson()
  expect_equal(dim(fitted(fit, summary = FALSE)), c(80, 118))
  # directly test the catordinal helper function
  mu <- fitted_catordinal(array(eta, dim = c(dim(eta), 3)), 
                          max_obs = 4, family = cumulative())
  expect_equal(dim(mu), c(80, 236, 4))
  # truncated helper functions
  # truncated continous models
  data <- list()
  mu <- fitted_trunc_gaussian(eta, lb = 0, ub = 10, x = fit, data = data)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_student(eta, lb = -Inf, ub = 15, x = fit, data = data)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_lognormal(eta, lb = 2, ub = 15, x = fit, data = data)
  expect_equal(dim(mu), c(80, 236))
  exp_eta <- exp(eta)
  mu <- fitted_trunc_gamma(exp_eta, lb = 1, ub = 7, x = fit)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_exponential(exp_eta, lb = 0, ub = Inf, x = fit)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_weibull(exp_eta, lb = -Inf, ub = Inf, x = fit)
  expect_equal(dim(mu), c(80, 236))
  # truncated discrete models
  data <- list(Y = sample(100, 10), trials = 1:10, N = 10)
  mu <- fitted_trunc_poisson(exp_eta, lb = 0, ub = 100, data = data)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_negbinomial(exp_eta, lb = 0, ub = 100, x = fit, 
                                 data = data)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_geometric(exp_eta, lb = 0, ub = 100, data = data)
  expect_equal(dim(mu), c(80, 236))
  mu <- fitted_trunc_binomial(ilink(eta, "logit"), lb = -Inf, ub = 100, 
                              data = data)
  expect_equal(dim(mu), c(80, 236))
})