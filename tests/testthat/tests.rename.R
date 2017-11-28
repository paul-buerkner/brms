test_that("make_index_names returns correct 1 and 2 dimensional indices", {
  expect_equal(make_index_names(rownames = 1:2), c("[1]", "[2]"))
  expect_equal(make_index_names(rownames = 1:2, colnames = 1:3, dim = 1), 
               c("[1]", "[2]"))
  expect_equal(make_index_names(rownames = c("a","b"), colnames = 1:3, dim = 2), 
               c("[a,1]", "[b,1]", "[a,2]", "[b,2]", "[a,3]", "[b,3]"))
})

test_that("rm_int_fe works as expected", {
  dat <- data.frame(y = 1:3,  y2 = 4:6, x = rnorm(3))
  code <- make_stancode(y ~ 1, data = dat)
  expect_equal(rm_int_fe("Intercept", code), character(0))
  code <- make_stancode(cbind(y, y2) ~ x, data = dat)
  expect_equal(rm_int_fe(c("Intercept", "x"), code, px = list(resp = "y")), "x")
  code <- make_stancode(y ~ x, data = dat, family = sratio())
  expect_equal(rm_int_fe(c("Intercept", "x"), code), "x")
  code <- make_stancode(y ~ 0 + intercept + x, data = dat)
  expect_equal(rm_int_fe(c("intercept", "x"), code), 
               c("intercept", "x"))
})
