test_that('set_inits produces the correct format', {
  res <- set_inits('normal(0, 1)', class = "Intercept", dpar = "mu")
  expect_s3_class(res, "data.frame")
  expect_s3_class(res, "brmsinits")
  res <- as.data.frame(res)
  expect_equal(res, data.frame(distribution = "normal(0, 1)",
                               class = "Intercept",
                               coef = "",
                               group = "",
                               dpar = "",
                               nlpar = ""))
})
