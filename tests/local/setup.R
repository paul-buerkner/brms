library(testthat)
library(brms)
set.seed(1234)

expect_ggplot <- function(object, ...) {
  testthat::expect_true(is(object, "ggplot"), ...)
}

expect_range <- function(object, lower = -Inf, upper = Inf, ...) {
  testthat::expect_true(all(object >= lower & object <= upper), ...)
}
