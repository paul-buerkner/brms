test_that("brm produces expected errors", {
  expect_error(brm(rating ~ treat + period + carry + cse(treat) + (1|subject), 
                   data = inhaler, family = cratio("logit")), 
               paste("Error occured for variables: treat"))
  expect_error(brm(rating ~ treat + period + carry + monotonic(carry),
                   data = inhaler, family = cratio("logit")), 
               paste("Error occured for variables: carry"))
})
