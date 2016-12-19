test_that("brm produces expected errors", {
  expect_error(brm(rating ~ treat + period + carry + cse(treat) + (1|subject), 
                   data = inhaler, family = cratio("logit")), 
               paste("Error occured for variables: treat"))
  expect_error(brm(rating ~ treat + period + carry + monotonic(carry),
                   data = inhaler, family = cratio("logit")), 
               paste("Error occured for variables: carry"))
})

test_that("check_brm_input returns correct warnings and errors", {
  x <- list(family = inverse.gaussian(), algorithm = "sampling")
  expect_warning(check_brm_input(x))
  x$family <- poisson("sqrt")
  expect_warning(check_brm_input(x))
})
