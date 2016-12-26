test_that("brmsformula validates formulas of non-linear parameters", {
  expect_error(bf(y ~ a, ~ 1, a ~ 1),
               "Additional formulas must be named")
  expect_error(bf(y ~ a^x, a.b ~ 1),
               "not contain dots or underscores")
  expect_error(bf(y ~ a^(x+b), a_b ~ 1),
               "not contain dots or underscores")
})

test_that("brmsformula validates formulas of auxiliary parameters", {
  expect_error(bf(y ~ a, ~ 1, sigma ~ 1),
               "Additional formulas must be named")
  expect_error(bf(y ~ a^x, a ~ 1),
               "The following parameter names are invalid: 'a'")
})

test_that("brmsformula does not change a 'brmsformula' object", {
  form <- bf(y ~ a, sigma ~ 1)
  expect_identical(form, bf(form))
  form <- bf(y ~ a, sigma ~ 1, a ~ x, nl = TRUE)
  expect_identical(form, bf(form))
})

test_that("brmsformula is backwards compatible", {
  expect_warning(form <- bf(y ~ a * exp(-b * x), 
                            nonlinear = a + b ~ 1),
                 "Argument 'nonlinear' is deprecated")
  expect_equivalent(pforms(form), list(a ~ 1, b ~ 1))
  expect_true(form[["nl"]])
  
  expect_warning(form <- bf(y ~ a * exp(-b * x), 
                            nonlinear = list(a ~ x, b ~ 1)),
                 "Argument 'nonlinear' is deprecated")
  expect_equivalent(pforms(form), list(a ~ x, b ~ 1))
  expect_true(form[["nl"]])
  
  form <- structure(y ~ x + z, sigma = sigma ~ x)
  class(form) <- c("brmsformula", "formula")
  form <- bf(form)
  expect_equal(form$formula, y ~ x + z)
  expect_equal(pforms(form), list(sigma = sigma ~ x))
  expect_true(!form[["nl"]])
  
  form <- structure(y ~ a * exp(-b * x),
                    nonlinear = list(a = a ~ x, b = b ~ 1))
  class(form) <- c("brmsformula", "formula")
  form <- bf(form)
  expect_equal(form$formula, y ~ a * exp(-b * x))
  expect_equal(pforms(form), list(a = a ~ x, b = b ~ 1))
  expect_true(form[["nl"]])
})
