test_that("brm produces expected errors", {
  expect_error(brm(rating~treat+period+carry+(1|subject), data = inhaler, 
                   partial = ~treat, family = cratio("logit")), 
              paste("Error occured for variables: treat"))
})