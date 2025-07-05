test_that("create_brm_call and call_only works", {
  # `call_only = T` creates a brm_call
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson() , call_only = T)
  # alternatively create_brm_call creates a brm_call without redefining parameters
  call_2  <- create_brm_call(count ~ zAge + zBase * Trt + (1|patient),
                             data = epilepsy, family = poisson())

  expect_equal(all.equal(call_1 , call_2), TRUE)
  c1 <- brm(call_1, call_only = TRUE)
  c2 <- brm(call_2, call_only = TRUE)
  expect_equal(all.equal(c1 , c2), TRUE)
})

test_that("brm_call objects are consistent", {
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T)
  call_2  <- create_brm_call(count ~ zBase * Trt + (1|patient),
                             data = epilepsy, family = poisson())
  expect_no_error(print(call_1))
  expect_no_error(print(call_2))
  expect_true(as.character(call_1$formula)[[3]] != as.character(call_2$formula)[[3]])
})

test_that("brm_call objects are consistent", {
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T)
  call_2 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson(), call_only = T)
  expect_no_error(print(call_1))
  expect_no_error(print(call_2))
  expect_true(as.character(call_1$formula)[[3]] == as.character(call_2$formula)[[3]])
})

test_that("brm_call objects are consistent", {
  # brm_call arguments will be checked and modify inside `brm()`function
  # therefore e.g. call_1$threads will ba a list containing threads etc.
  call_1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
                data = epilepsy, family = poisson() , call_only = T)
  expect_no_error(print(call_1))
  expect_true(all(is.null(call_1$data2) , call_1$iter == 2000, call_1$normalize,
                   call_1$warmup == 1000, is.null(call_1$threads$threads)))
})
