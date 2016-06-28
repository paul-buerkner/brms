test_that("nonlinear_predictor runs without errors", {
  # use brmsfit_example as the basis for constructing
  # a simple model that can be used to test nonlinear_predictor
  # more complicated non-linear models are only tested locally
  # to keep the installation time of brms short
  fit <- rename_pars(brmsfit_example)
  fit$formula <- bf(count ~ alpha - exp(Trt / beta),
                    nonlinear = list(alpha ~ 1, beta ~ 1))
  fit <- add_samples(fit, "b_alpha_Intercept")
  fit <- add_samples(fit, "b_beta_Intercept")
  draws <- extract_draws(fit)
  expect_equal(dim(nonlinear_predictor(draws)), 
               c(Nsamples(fit), nobs(fit)))
  newdata <- data.frame(count = rpois(3, 10), visit = 1:3, patient = 10,
                        Trt = rnorm(3), Age = 0, Exp = 1)
  draws <- extract_draws(fit, newdata = newdata)
  expect_equal(dim(nonlinear_predictor(draws)), 
               c(Nsamples(fit), nrow(newdata)))
})