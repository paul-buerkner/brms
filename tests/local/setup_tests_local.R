library(testthat)
library(brms)
options(mc.cores = 2)
ggplot2::theme_set(theme_default())
# prevents rstan related crashes
options(brms.backend = "cmdstanr")
set.seed(1234)

expect_ggplot <- function(object, ...) {
  testthat::expect_true(is(object, "ggplot"), ...)
}

expect_range <- function(object, lower = -Inf, upper = Inf, ...) {
  testthat::expect_true(all(object >= lower & object <= upper), ...)
}

SW <- suppressWarnings
SM <- suppressMessages

context("local tests")
