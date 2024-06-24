context("Tests for priorsense support")

skip_on_cran()

require(priorsense)

fit1 <- rename_pars(brms:::brmsfit_example1)

test_that("create_priorsense_data returns expected output structure", {
  psd1 <- create_priorsense_data(fit1)
  expect_s3_class(psd1$draws, "draws")
  expect_s3_class(psd1$fit, "brmsfit")
  expect_s3_class(psd1$log_lik, "draws")
  expect_s3_class(psd1$log_prior, "draws")
  expect_true(is.function(psd1$log_lik_fn))
  expect_true(is.function(psd1$log_prior_fn))
  expect_true(is.function(psd1$log_ratio_fn))
})

test_that("powerscale returns without error", {
  expect_no_error(powerscale(fit1, component = "prior", alpha = 0.8))
  expect_no_error(powerscale(fit1, component = "likelihood", alpha = 1.1))
})

test_that("powerscale_sensitivity returns without error", {
  expect_no_error(powerscale_sensitivity(fit1))
})
