test_that("brm hash method works", {

  # same all
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson() )

  expect_true( c1 == c2 )

  # different formula
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~   zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson() )

  expect_false( c1 == c2 )

  # different data
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~   zBase * Trt + (1|patient),
                                data = epilepsy[-c(1), ], family = poisson() )

  expect_false( c1 == c2 )



  a <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                              data = epilepsy[-c(1) , ], family = poisson() )

  b <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                              data = epilepsy[-c(1) , ], family = poisson() )

  expect_true( a == b )

  a <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                              data = epilepsy[-c(2) , ], family = poisson() )

  b <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                              data = epilepsy[-c(1) , ], family = poisson() )

  expect_false( a == b )


})

# tests/testthat/test-hash_brm_call_master.R
# -------------------------------------------------------------------
#  Tests for the top-level hash_brm_call_master() wrapper
# -------------------------------------------------------------------

data(epilepsy, package = "brms")  # dataset used in several tests

make_args <- function(fam = poisson()) {
  list(
    formula = count ~ zAge + zBase,
    data    = epilepsy,
    family  = fam
  )
}

test_that("returns a non-empty character string", {
  key <- hash_brm_call_master(make_args())
  expect_true(is.character(key) && length(key) == 1L && nchar(key) > 0)
})

test_that("order of names in args_list does not change hash", {
  key1 <- hash_brm_call_master(make_args())

  # shuffle the order of elements
  shuffled <- make_args()[c("family", "data", "formula")]
  key2 <- hash_brm_call_master(shuffled)

  expect_identical(key1, key2)
})

test_that("changing a model-defining argument changes the hash", {
  key_pois <- hash_brm_call_master(make_args(poisson()))
  key_gaus <- hash_brm_call_master(make_args(gaussian()))
  expect_false(key_pois == key_gaus)
})

test_that("unnamed or partially named lists trigger error", {
  bad1 <- list(1, 2, 3)                     # no names
  bad2 <- list(formula = count ~ zAge, 2)   # one missing name
  expect_error(hash_brm_call_master(bad1), "named")
  expect_error(hash_brm_call_master(bad2), "named")
})

test_that("hash is stable across repeated identical calls", {
  k1 <- hash_brm_call_master(make_args())
  k2 <- hash_brm_call_master(make_args())
  expect_identical(k1, k2)
})

test_that("empty list gives informative error", {
  expect_error(hash_brm_call_master(list()), "named")
})

