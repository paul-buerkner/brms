test_that("plot doesn't throw errors", {
  fit <- rename_pars(brmsfit_example)
  expect_silent(p <- plot(fit, do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^b", do_plot = FALSE))
  expect_silent(p <- plot(fit, pars = "^sd", do_plot = FALSE))
  expect_error(p <- plot(fit, pars = "123", do_plot = FALSE),
               "No valid parameters selected")
})

test_that("stanplot and pairs doesn't throw errors", {
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
  #expect_silent(p <- stanplot(fit, type = "pairs", quiet = TRUE,
  #                            pars = parnames(fit)[2:3]))
  #expect_silent(p <- stanplot(fit, type = "diag", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "rhat", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ess", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "mcse", quiet = TRUE))
  expect_silent(p <- stanplot(fit, type = "ac", quiet = TRUE))
  # tests for pairs (deactivated for the moment)
  #expect_silent(p <- pairs(fit, pars = parnames(fit)[1:3]))
})