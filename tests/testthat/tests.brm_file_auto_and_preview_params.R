
test_that("preview parameter for brm function works", {

  # This test will be faster since it will use preview parameter and get a list
  # evaluated parameters for the call and the hash value for the call.

  fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                      data = epilepsy, family = poisson() , file_auto = TRUE , preview = TRUE  )

  fit1_gaussian <- brm(count ~ zAge + zBase * Trt + (1|patient),
              data = epilepsy, family = gaussian() , file_auto = TRUE , preview = TRUE  )

  # override any parameter of a brm_call_preview with realize
  fit1_gaussian2 <- realize( fit1 , family= gaussian() , preview= TRUE )

  expect_equal(  fit1_gaussian$hash , fit1_gaussian2$hash )

  fit2 <- brm(count ~ zAge + zBase * Trt + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE , preview = TRUE)

  fit3 <- brm(count ~ zAge + zBase  + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE  , preview = TRUE )

  fit4 <- brm(count ~ zAge + zBase  + (1|patient),
              data = epilepsy, family = poisson() , file_auto = FALSE , preview = TRUE )





  # fit2 will return same result with fit1 cause file_auto is TRUE
  expect_equal(  fit1$hash  ,  fit2$hash )
  # fit1 and fit3 must have different hashes
  expect_false(  fit1$hash  ==  fit3$hash )

  # fit3 and fit4 has same parameters that has an effect on the result
  expect_equal(  fit3$hash , fit4$hash )

  # check if family parameter also changes the hash
  expect_false(  fit1$hash  ==  fit1_gaussian$hash )

})

test_that("realize function evaluates brm_call_preview", {

  #  when preview parameter is TRUE brm function will return a brm_call_preview
  #  user will be able to get real result with realize function as below
  old_val <- getOption("brms.cache_folder", default= NULL)
  cache_folder <- "ignore_caches"
  options(brms.cache_folder = cache_folder  )
  if(!dir.exists(cache_folder))
    dir.create(cache_folder )
  on.exit({
    if (is.null(old_val)) {
      options(brms.cache_folder = NULL)
    } else {
      options(brms.cache_folder = old_val)
    }
  }, add = TRUE)


  fit1_preview <- brm(count ~ zAge + zBase * Trt + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE   , preview = TRUE  )

  fit2_preview <- brm(count ~ zAge + zBase  + (1|patient),
              data = epilepsy, family = poisson() , file_auto = TRUE , preview = TRUE)

  expect_false(  fit1_preview$hash  ==   fit2_preview$hash )

  fit1_real <-  realize(fit1_preview)
  fit2_real <-  realize(fit2_preview)
  expect_equal(  fit1_real$run_info$hash  ,  fit1_real$run_info$hash )
})


test_that("file_auto option works", {

  old_val <- getOption("brms.cache_folder", default= NULL)
  cache_folder <- "ignore_caches"
  options(brms.cache_folder = cache_folder  )
  if(!dir.exists(cache_folder))
    dir.create(cache_folder )
  on.exit({
    if (is.null(old_val)) {
      options(brms.cache_folder = NULL)
    } else {
      options(brms.cache_folder = old_val)
    }
  }, add = TRUE)



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

