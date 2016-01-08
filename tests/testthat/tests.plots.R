test_that("plot doesn't throw errors", {
  fit <- rename_pars(brmsfit_example)
  expect_silent(p <- plot(fit, do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^b", do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^sd", do_plot = FALSE))
  expect_error(p <- plot(fit, pars = "123", do_plot = FALSE),
               "No valid parameters selected")
})

test_that("stanplot and pairs works correctly", {
  fit <- rename_pars(brmsfit_example)
  # tests for stanplot
  expect_silent(p <- stanplot(fit, quiet = TRUE))
  expect_silent(p <- stanplot(fit, pars = "^b", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "trace", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "hist", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "dens", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "scat", quiet = TRUE,
                              pars = parnames(fit)[2:3], 
                              exact_match = TRUE))
  #expect_silent(p <- stanplot(fit, type = "diag", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "rhat", pars = "^b_",
                              quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ess", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "mcse", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ac", quiet = TRUE))
  expect_identical(pairs(fit, pars = parnames(fit)[1:3]), NULL)
  # warning occurs somewhere in rstan
  expect_warning(stanplot(fit, type = "par", pars = "^b_Intercept$"))
  expect_warning(p <- stanplot(fit, type = "par", pars = "^b_"),
                 "stan_par expects a single parameter name")
  expect_error(stanplot(fit, type = "density"), "Invalid plot type")
})