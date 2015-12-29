test_that("brm produces expected errors", {
  expect_error(brm(rating~treat+period+carry+(1|subject), data = inhaler, 
                   partial = ~treat, family = c("cratio", "logit")), 
              paste("Variables cannot be modeled as fixed", 
                    "and partial effects at the same time.", 
                    "Error occured for variables: treat"))
})