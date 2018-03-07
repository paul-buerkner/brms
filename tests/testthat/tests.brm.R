# calling context() avoids a strange bug in testthat 2.0.0
# cannot actually run brms models in tests as it takes way too long
context("Tests for brms error messages")

test_that("brm produces expected errors", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
  
  # formula parsing
  expect_error(brm(bf(y ~ a^x, a + b ~ 1, nl = TRUE), dat),
               "missing in formula: 'b'")
  expect_error(brm(~ x + (1|g), dat), 
               "response variable is missing")
  expect_error(brm(bf(y ~ a, nl = TRUE)),
               "No non-linear parameters specified")
  expect_error(brm(bf(y ~ a, a ~ 1, nl = TRUE), family = acat()),
               "Non-linear formulas are not yet allowed for this family")
  expect_error(brm(bf(y ~ 0 + a), dat, family = cumulative()),
               "Cannot remove the intercept in an ordinal model")
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
  
  # autocorrelation
  expect_error(brm(y ~ 1, dat, autocor = cor_ar(~x+y|g)), 
               "Autocorrelation structures may only contain 1 time variable")
  expect_error(brm(y ~ 1, dat, autocor = cor_ar(x~t1|g1)), 
               "Autocorrelation formulas must be one-sided")
  expect_error(brm(y ~ 1, dat, autocor = cor_ar(~1|g1/g2)), 
               paste("Illegal grouping term: g1/g2"))
  expect_error(brm(y ~ 1, dat, poisson(), autocor = cor_ma(~x)),
               "not implemented for family 'poisson'")
  
  # ordinal models
  expect_error(brm(rating ~ treat + (cs(period)|subject),
                   data = inhaler, family = categorical()), 
               "Category specific effects are only meaningful")
  
  # families and links
  expect_error(brm(y ~ x, dat, family = gaussian("logit")), 
               "'logit' is not a supported link for family 'gaussian'")
  expect_error(brm(y ~ x, dat, family = poisson("inverse")), 
               "'inverse' is not a supported link for family 'poisson'")
  expect_error(brm(y ~ x, dat, family = c("weibull", "sqrt")),
               "'sqrt' is not a supported link for family 'weibull'")
  expect_error(brm(y ~ x, dat, family = c("categorical", "probit")),
               "'probit' is not a supported link for family 'categorical'")
  expect_error(brm(y ~ x, dat, family = "ordinal"),
              "ordinal is not a supported family")
})
