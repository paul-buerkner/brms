# calling context() avoids a strange bug in testthat 2.0.0
# cannot actually run brms models in tests as it takes way too long
context("Tests for brms error messages")

test_that("brm works fully with mock backend", {
  skip_on_cran()
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  # Positive control - forced error gets thrown and propagated
  expect_error(brm(y ~ x + (1|g), dat, backend = "mock",
                   stan_model_args = list(compile_error = "Test error")),
               "Test error")

  # Positive control - bad Stan code from stanvars gets an error
  # test passes but prints output that is somehow impossible to suppress
  # expect_error(suppressMessages(
  #    brm(y ~ x + (1|g), dat, backend = "mock",
  #        stanvars = stanvar(scode = "invalid;", block = "model"))
  # ))

  # Testing some models
  mock_fit <- brm(y ~ x + (1|g), dat, mock_fit = 1, backend = "mock", rename = FALSE)
  expect_equal(mock_fit$fit, 1)
})

test_that("brm(file = xx) works fully with mock backend", {
  skip_on_cran()

  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  file <- tempfile(fileext = ".rds")
  mock_fit1 <- brm(y ~ x + (1|g), dat, mock_fit = "stored", backend = "mock",
                   rename = FALSE, file = file)
  expect_true(file.exists(file))

  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")

  # In default settings, even using different data/model should result in the
  # model being loaded from file
  changed_data <- dat[1:8, ]
  mock_fit2 <- brm(y ~ x + 0, changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")

  # Now test using file_refit = "on_change" which should be more clever
  # No change
  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file)
  expect_equal(mock_fit2$fit, "stored")


  # Change data, but not code
  mock_fit2 <- brm(y ~ x + (1|g), changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change")
  expect_equal(mock_fit2$fit, "new")

  # Change code but not data
  mock_fit2 <- brm(y ~ x + (1|g), dat, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change",
                   prior = prior(normal(0,2), class = sd))
  expect_equal(mock_fit2$fit, "new")

  # Change both
  mock_fit2 <- brm(y ~ x + 0, changed_data, mock_fit = "new", backend = "mock",
                   rename = FALSE, file = file, file_refit = "on_change")
  expect_equal(mock_fit2$fit, "new")
})


test_that("brm produces expected errors", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))

  # formula parsing
  expect_error(brm(~ x + (1|g), dat, file = "test"),
               "Response variable is missing")
  expect_error(brm(bf(y ~ a, nl = TRUE)),
               "No non-linear parameters specified")
  expect_error(brm(bf(y | se(sei) ~ x, sigma ~ x), dat),
               "Cannot predict or fix 'sigma' in this model")
  expect_error(brm(y | se(sei) ~ x, dat, family = weibull()),
               "Argument 'se' is not supported for family")
  expect_error(brm(y | se(sei) + se(sei2) ~ x, dat, family = gaussian()),
               "Each addition argument may only be defined once")
  expect_error(brm(y | abc(sei) ~ x, family = gaussian()),
               "The following addition terms are invalid:\n'abc(sei)'",
               fixed = TRUE)
  expect_error(brm(y | disp(sei) ~ x, dat, family = gaussian()),
               "The following addition terms are invalid:")
  expect_error(brm(bf(y ~ x, shape ~ x), family = gaussian()),
               "The parameter 'shape' is not a valid distributional")
  expect_error(brm(y ~ x + (1|abc|g/x), dat),
               "Can only combine group-level terms")
  expect_error(brm(y ~ x + (1|g) + (x|g), dat),
               "Duplicated group-level effects are not allowed")
  expect_error(brm(y~mo(g)*t2(x), dat), fixed = TRUE,
               "The term 'mo(g):t2(x)' is invalid")
  expect_error(brm(y~x*cs(g), dat), fixed = TRUE,
               "The term 'x:cs(g)' is invalid")
  expect_error(brm(y~me(x, 2 * g)*me(x, g), dat),
               "Variable 'x' is used in different calls to 'me'")
  expect_error(brm(y ~ 1 + set_rescor(TRUE), data = dat),
               "Function 'set_rescor' should not be part")

  # autocorrelation
  expect_error(brm(y ~ ar(x+y, g), dat),
               "Cannot coerce 'x \\+ y' to a single variable name")
  expect_error(brm(y ~ ar(gr = g1/g2), dat),
               "Illegal grouping term 'g1/g2'")
  expect_error(brm(y ~ ma(x), dat, poisson()),
               "Please set cov = TRUE")
  expect_error(brm(bf(y ~ 1) + arma(x), dat),
               "Autocorrelation terms can only be specified")

  # ordinal models
  expect_error(brm(rating ~ treat + (cs(period)|subject),
                   data = inhaler, family = categorical()),
               "Category specific effects are not supported")

  # families and links
  expect_error(brm(y ~ x, dat, family = poisson("inverse")),
               "'inverse' is not a supported link for family 'poisson'")
  expect_error(brm(y ~ x, dat, family = c("weibull", "sqrt")),
               "'sqrt' is not a supported link for family 'weibull'")
  expect_error(brm(y ~ x, dat, family = c("categorical", "probit")),
               "'probit' is not a supported link for family 'categorical'")
  expect_error(brm(y ~ x, dat, family = "ordinal"),
              "ordinal is not a supported family")
})

