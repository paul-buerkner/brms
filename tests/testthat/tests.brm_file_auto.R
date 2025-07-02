
test_that("hash function check" , {
  aa1 =   hash_func_for_testing( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  aa2 =   hash_func_for_testing( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  expect_equal( aa1 , aa2 )

})



# ─────────────────────────────────────────────────────────────────────────
test_that("List argument order is normalised", {
  f <- bf(mpg ~ wt)

  h1 <- hash_func_for_testing(formula = f, data = mtcars, family = gaussian(), iter = 2000)
  h2 <- hash_func_for_testing(iter = 2000, family = gaussian(), data = mtcars, formula = f)

  expect_identical(h1, h2)
})


