# tests/testthat/test-stop2.R

test_that("stop2() throws a brms_error with simple message", {
  expect_error(stop2("Something went wrong"), class = "brms_error")
})

test_that("stop2() throws a brms_error with simple message", {
  value <- 42
  expect_error(stop2("Unexpected value {value}"), class = "brms_error")
})

test_that("stop2() supports custom subclasses", {
  err <- expect_error(
    stop2("Some issue", .subclass = "brms_custom"),
    class = "brms_custom"
  )
  expect_true(inherits(err, "brms_error"))
  expect_true(inherits(err, "brms_custom"))
})

test_that("stop2() supports custom subclasses", {
  algorithm = "SomeAlgorithm"
  err <- expect_error(
    stop2("Algorithm '", algorithm, "' is not supported." ,.subclass = "brms_algorithm" ),
    class = "brms_algorithm"
  )
  expect_true(inherits(err, "brms_algorithm"))
})

test_that("stop2() respects an explicit call argument", {
   err <- expect_error(
    stop2("Cannot remove the intercept in an ordinal model." , .subclass = "brms_ordinal_model"),
    class = "brms_ordinal_model"
  )
  expect_true(inherits(err, "brms_ordinal_model"))
})

test_that("stop2() respects an explicit call argument", {
  fake_call <- quote(fake_func("a", 123))
  err <- expect_error(stop2("error occurred", call = fake_call))
  expect_equal(conditionCall(err), fake_call)
})
