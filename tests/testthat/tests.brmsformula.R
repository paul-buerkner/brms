context("Tests for brmsformula")

test_that("brmsformula validates formulas of non-linear parameters", {
  expect_error(
    bf(y ~ a, ~1, a ~ 1),
    "Additional formulas must be named"
  )
  expect_error(
    bf(y ~ a^x, a.b ~ 1),
    "not contain dots or underscores"
  )
  expect_error(
    bf(y ~ a^(x + b), a_b ~ 1),
    "not contain dots or underscores"
  )
})

test_that("brmsformula validates formulas of auxiliary parameters", {
  expect_error(
    bf(y ~ a, ~1, sigma ~ 1),
    "Additional formulas must be named"
  )
})

test_that("brmsformula detects use if '~~'", {
  # checks fix of issue #749
  expect_error(bf(y ~ ~x), "~~")
})

test_that("brmsformula does not change a 'brmsformula' object", {
  form <- bf(y ~ a, sigma ~ 1)
  expect_identical(form, bf(form))
  form <- bf(y ~ a, sigma ~ 1, a ~ x, nl = TRUE)
  expect_identical(form, bf(form))
})

test_that("brmsformula detects auxiliary parameter equations", {
  expect_error(
    bf(y ~ x, sigma1 = "sigmaa2"),
    "Can only equate parameters of the same class"
  )
  expect_error(
    bf(y ~ x, mu3 = "mu2"),
    "Equating parameters of class 'mu' is not allowed"
  )
  expect_error(
    bf(y ~ x, sigma1 = "sigma1"),
    "Equating 'sigma1' with itself is not meaningful"
  )
  expect_error(
    bf(y ~ x, shape1 ~ x, shape2 = "shape1"),
    "Cannot use predicted parameters on the right-hand side"
  )
  expect_error(
    bf(y ~ x, shape1 = "shape3", shape2 = "shape1"),
    "Cannot use fixed parameters on the right-hand side"
  )
})

test_that("update_adterms works correctly", {
  form <- y | trials(size) ~ x
  expect_equal(
    update_adterms(form, ~ trials(10)),
    y | trials(10) ~ x
  )
  expect_equal(
    update_adterms(form, ~ weights(w)),
    y | trials(size) + weights(w) ~ x
  )
  expect_equal(
    update_adterms(form, ~ weights(w), action = "replace"),
    y | weights(w) ~ x
  )
  expect_equal(
    update_adterms(y ~ x, ~ trials(10)),
    y | trials(10) ~ x
  )
})
