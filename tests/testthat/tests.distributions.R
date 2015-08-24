test_that("Test that student distribution functions work correctly", {
  expect_equal(integrate(dstudent, -100, 100, df = 15, mu = 10, sigma = 5)$value, 1)
  expect_equal(dstudent(1, df = 10, mu = 0, sigma = 5), dt(1/5, df = 10)/5)
  expect_equal(pstudent(2, df = 20, mu = 2, sigma = 0.4), pt(0, df = 20))
  expect_equal(qstudent(0.7, df = 5, mu = 2, sigma = 3), 2 + 3*qt(0.7, df = 5))
})