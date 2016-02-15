test_that("nonlinear_predictor runs without errors", {
  # use brmsfit_example as the basis for constructing
  # a simple model that can be used to test nonlinear_predictor
  # more complicated non-linear models are only tested locally
  # to keep the installation time of brms short
  fit <- rename_pars(brmsfit_example)
  fit$formula <- count ~ alpha - exp(Trt_c / beta)
  fit$nonlinear <- list(alpha ~ 1, beta ~ 1)
  fit <- add_samples(fit, "b_alpha_Intercept")
  fit <- add_samples(fit, "b_beta_Intercept")
  C <- model.frame(fit)[, "Trt_c", drop = FALSE]
  expect_equal(dim(nonlinear_predictor(fit, C = C)), 
               c(Nsamples(fit), nobs(fit)))
  newdata <- data.frame(count = rpois(3, 10), visit = 1:3, patient = 10)
  C <- structure(matrix(rnorm(3)), dimnames = list(NULL, "Trt_c"))
  expect_equal(dim(nonlinear_predictor(fit, C = C, newdata = newdata)), 
               c(Nsamples(fit), nrow(newdata)))
})