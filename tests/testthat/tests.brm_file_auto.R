test_that("file_auto option works", {
  fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE )
  fit2 <- brm(count ~ zAge + zBase * Trt + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE )
  fit3 <- brm(count ~ zAge + zBase  + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE   )
  fit4 <- brm(count ~ zAge + zBase  + (1|patient),
              data = epilepsy, family = poisson() , file_auto = FALSE  )
  fit1_gaussian <- brm(count ~ zAge + zBase * Trt + (1|patient),
                       data = epilepsy, family = gaussian() , file_auto = TRUE )

  # fit2 will return same result with fit1 cause file_auto is TRUE
  expect_equal(  fit1$run_info$hash  ,  fit2$run_info$hash )
  # fit1 and fit3 must have different hashes
  expect_false(  fit1$run_info$hash  ==  fit3$run_info$hash )
  # fit3 and fit4 has same parameters that has an effect on the result
  expect_equal(  fit3$run_info$hash , fit4$run_info$hash )
  # but fit4 should trigger a refit and should have a different start time
  expect_false(  fit3$run_info$start  ==  fit4$run_info$start )
  # check if family parameter also changes the hash
  expect_false(  fit1$run_info$hash  ==  fit1_gaussian$run_info$hash )
})