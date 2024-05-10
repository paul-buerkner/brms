source("tests/local/setup_tests_local.R")

data <- epilepsy
data$cat <- factor(sample(1:5, nrow(data), replace = TRUE))

test_that('gaussian model runs with set_inits', {
  formula <- bf(count ~ cat + Age,
                sigma ~ cat + Age)

  inits <- set_inits('normal(0, 1)', class = "Intercept", dpar = "mu") +
    set_inits('uniform(-0.1, 0.1)', class = "b", dpar = "sigma")

  out <- capture_messages(fit <- brm(formula, data = data, init = inits, refresh = 0, chains = 2))
  out <- paste0(out, collapse = "\n")
  expect_true(grepl("Missing init values for the following parameters:", out))
  expect_true(grepl(": b, Intercept_sigma\n", out))
  fit_init <- fit$stan_args$init
  expect_length(fit_init, 2)
  expect_equal(names(fit_init[[1]]), c("Intercept", "b_sigma"))
  expect_range(fit_init[[1]]$b_sigma, -0.1, 0.1)
})


test_that('poisson model runs with set_inits', {
  formula <- bf(count ~ cat + Age)

  inits <- set_inits('normal(0, 1)', class = "Intercept", dpar = "mu") +
    set_inits('uniform(-0.1, 0.1)', class = "b", dpar = "mu")

  out <- capture_messages(fit <- brm(formula, data = data, family = poisson(),
                                     init = inits, refresh = 0, chains = 2))
  out <- paste0(out, collapse = "\n")
  expect_false(grepl("Missing init values for the following parameters:", out))
  fit_init <- fit$stan_args$init
  expect_length(fit_init, 2)
  expect_equal(names(fit_init[[1]]), c("Intercept", "b"))
  expect_range(fit_init[[1]]$b, -0.1, 0.1)
})
