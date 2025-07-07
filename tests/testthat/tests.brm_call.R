test_that("create_brm_call and call_only works", {

  file_auto <- TRUE
  # `call_only = T` creates a brm_call
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T, file_auto = file_auto)
  # alternatively create_brm_call creates a brm_call without redefining parameters
  call_2  <- create_brm_call(count ~ zAge + zBase * Trt + (1|patient),
                             data = epilepsy, family = poisson(),
                             file_auto = file_auto)

  expect_equal(all.equal(call_1 , call_2), TRUE)
  c1 <- brm(call_1, call_only = TRUE)
  c2 <- brm(call_2, call_only = TRUE)
  expect_equal(all.equal(c1 , c2), TRUE)
})

test_that("brm_call objects are consistent", {
  file_auto <- TRUE
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T, file_auto = file_auto)
  call_2  <- create_brm_call(count ~ zBase * Trt + (1|patient),
                             data = epilepsy, family = poisson(), file_auto = file_auto)
  expect_no_error(print(call_1))
  expect_no_error(print(call_2))
  expect_true(as.character(call_1$formula)[[3]] != as.character(call_2$formula)[[3]])
})

test_that("brm_call objects are consistent", {
  file_auto <- TRUE
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T, file_auto = file_auto)
  call_2 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T, file_auto = file_auto)
  expect_no_error(print(call_1))
  expect_no_error(print(call_2))
  expect_true(as.character(call_1$formula)[[3]] == as.character(call_2$formula)[[3]])
})

