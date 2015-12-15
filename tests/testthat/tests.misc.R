test_that("Test that is.formula is TRUE for formulas and otherwise FALSE", {
  expect_equal(is.formula(y~1), TRUE)
  expect_equal(is.formula("a"), FALSE)
  expect_equal(is.formula(list(y~1, ~1)), TRUE)
  expect_equal(is.formula(list(y~1,1)), TRUE)
  expect_equal(is.formula(list("a",1)), FALSE)
  expect_equal(is.formula(list(y~1, ~1), or = FALSE), TRUE)
  expect_equal(is.formula(list(y~1,1), or = FALSE), FALSE)
})

test_that("Test that rmNULL removes all NULL entries", {
  expect_equal(rmNULL(list(a = NULL, b = 1, c = list(NULL, 1))),
               list(b = 1, c = list(1)))
  expect_equal(rmNULL(list(a = NULL, b = 1, c = NULL)),
               list(b = 1))
})

test_that("Test that rmNum remove all numeric entries", {
  expect_equal(rmNum(list(1, "a", 2.3, "b")), list("a","b"))
  expect_equal(rmNum(list(x = 1.5, y = "abc", z = pi)), list(y = "abc"))
})

test_that("Test that forumla2string performs correct conversion", {
  expect_error(formula2string("y~x"))
  expect_equal(formula2string(y ~ x + c), "y~x+c")
  expect_equal(formula2string(abc ~ x + cd, rm = c(3,2)), "~x+")
})

test_that("Test that collapse_lists performs correct collapsing after names", {
  x <- list(a = "a <- ", b = "b <- ")
  y <- list(b = "cauchy(1,2)", c = "normal(0,1)", a = "gamma(1,1)")
  expect_equal(collapse_lists(list()), list())
  expect_equal(collapse_lists(list(x, y)), 
               list(a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)", c = "normal(0,1)"))
  expect_equal(collapse_lists(list(c(x, c = "c <- "), y)),
               list(a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)", c = "c <- normal(0,1)"))
})

test_that("Test that nlist works correctly", {
  x <- 1
  y <- 2:3
  exlist <- list(x = x, y = y)
  expect_equal(nlist(x = x, y = y), exlist)
  expect_equal(nlist(x, y), exlist)
  expect_equal(nlist(x = x, y), exlist)
})

test_that("Test that convenience functions for model families work correctly", {
  expect_true(use_real(gaussian()))
  expect_true(!use_real("poisson"))
  expect_true(use_int(binomial()))
  expect_true(has_trials("zero_inflated_binomial"))
  expect_true(has_cat("acat"))
  expect_true(has_sigma(student()))
  expect_true(!has_sigma("cauchy", se = TRUE))
  expect_true(has_sigma("cauchy", se = TRUE, autocor = cor_ar()))
})

test_that("Test that check_intercept updates FE names", {
  expect_equal(check_intercept(c("Intercept", "x", "z")),
               list(names = c("x", "z"), has_intercept = TRUE))
  expect_equal(check_intercept(c("x", "z")),
               list(names = c("x", "z"), has_intercept = FALSE))
})