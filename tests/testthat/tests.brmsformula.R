test_that("brmsformula validates formulas of non-linear parameters", {
  expect_error(bf(y ~ a, ~ 1, a ~ 1),
               "Additional formulas must be named")
  expect_error(bf(y ~ a^x, a.b ~ 1),
               "not contain dots or underscores")
  expect_error(bf(y ~ a^(x+b), a_b ~ 1),
               "not contain dots or underscores")
})

test_that("brmsformula is backwards compatible", {
  expect_warning(form <- bf(y ~ a * exp(-b * x), 
                            nonlinear = a + b ~ 1),
                 "Argument 'nonlinear' is deprecated")
  expect_equivalent(pforms(form), list(a ~ 1, b ~ 1))
  expect_warning(form <- bf(y ~ a * exp(-b * x), 
                            nonlinear = list(a ~ x, b ~ 1)),
                 "Argument 'nonlinear' is deprecated")
  expect_equivalent(pforms(form), list(a ~ x, b ~ 1))
})

test_that("brmsformula validates formulas of auxiliary parameters", {
  expect_error(bf(y ~ a, ~ 1, sigma ~ 1),
               "Additional formulas must be named")
  expect_error(bf(y ~ a^x, a ~ 1),
               "The following parameter names are invalid: 'a'")
})
