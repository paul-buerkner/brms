test_that("brmsformula validates formulas of non-linear parameters", {
  expect_error(bf(y ~ a, nonlinear = list(~ 1, a ~ 1)),
               "Additional formulas must be named")
  expect_error(bf(y ~ exp(-x), nonlinear = list(a + b ~ 1)),
               "LHS of additional formulas must contain exactly one variable")
  expect_error(bf(y ~ a^x, nonlinear = list(a.b ~ 1)),
               "not contain dots or underscores")
  expect_error(bf(y ~ a^(x+b), nonlinear = a_b ~ 1),
               "not contain dots or underscores")
})

test_that("brmsformula validates formulas of auxiliary parameters", {
  expect_error(bf(y ~ a, ~ 1, sigma ~ 1),
               "Additional formulas must be named")
  expect_error(bf(y ~ exp(-x), a + b ~ 1),
               "LHS of additional formulas must contain exactly one variable")
  expect_error(bf(y ~ a^x, a ~ 1),
               "The following argument names were invalid: a")
})

test_that("nonlinear2list works correctly", {
  expect_equal(nonlinear2list(a ~ 1), list(a = a ~ 1))
  expect_equal(nonlinear2list(a + alpha ~ x + (x|g)), 
               list(a = a ~ x + (x|g), alpha = alpha ~ x + (x|g)))
  expect_equal(nonlinear2list(list(a ~ 1, b ~ 1 + z)),
               list(a = a ~ 1, b = b ~ 1 + z))
  expect_equal(nonlinear2list(NULL), NULL)
  expect_error(nonlinear2list(1), "invalid formula")
})