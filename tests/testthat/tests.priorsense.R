context("Tests for priorsense support")

SW <- suppressWarnings

skip_on_cran()

require(priorsense)

fit1 <- rename_pars(brms:::brmsfit_example1)

fit2 <- rename_pars(brms:::brmsfit_example7)

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
  expect_no_error(SW(powerscale(fit1, component = "prior", alpha = 0.8)))
  expect_no_error(SW(powerscale(fit1, component = "likelihood", alpha = 1.1)))
})

test_that("powerscale_sensitivity returns without error", {
  expect_no_error(SW(powerscale_sensitivity(fit1)))
})

test_that("prior tags correctly go to priorsense data", {

  psd <- create_priorsense_data.brmsfit(fit2)

  expect_setequal(
    variables(psd$log_prior),
    c(
      "lprior",
      "lprior_prior_tag1",
      "lprior_prior_tag2",
      "lprior_prior_tag3",
      "lprior_prior_tag4",
      "lprior_prior_tag5"
    )
  )
})
