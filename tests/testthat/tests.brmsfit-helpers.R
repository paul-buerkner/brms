test_that("Test that first_greater returns expected results", {
  set.seed(12345)
  A <- cbind(1:10, 11:20, 21:30)
  x <- sample(1:25, size = 10, replace = TRUE)
  expect_equal(first_greater(A, x), c(3,3,3,3,2,1,2,2,2,3))
  expect_equal(first_greater(A, x, i = 2), c(3,3,3,3,2,2,2,2,2,3))
})

test_that("Test that array2list performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(matrix(1:9,3,3), matrix(10:18,3,3), matrix(19:27,3,3))
  expect_equal(array2list(A), B)
})

test_that("Test that list2array performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(array(1:9, dim = c(3,3,1)), array(10:18, dim = c(3,3,1)), array(19:27, dim = c(3,3,1)))
  C <- list(array(1:9, dim = c(3,3)), array(10:18, dim = c(3,3)), array(19:27, dim = c(3,3)))
  expect_equal(list2array(B), A)
  expect_equal(list2array(C), A)
})

test_that("Test that ilink(x, probit) and ilink(x, probit_approx) produce similar results", {
  all.equal(ilink(-10:10, "probit"), ilink(-10:10, "probit_approx"), tolerance = 1e-3)
})

test_that("Test that get_cornames returns desired correlation names", {
  names <- c("Intercept", "x", "y")
  expect_equal(get_cornames(names), c("cor(Intercept,x)", "cor(Intercept,y)", "cor(x,y)"))
  expect_equal(get_cornames(names, brackets = FALSE), 
               c("cor_Intercept_x", "cor_Intercept_y", "cor_x_y"))
  expect_equal(get_cornames(names, type = "rescor"),
               c("rescor(Intercept,x)", "rescor(Intercept,y)", "rescor(x,y)"))
  expect_equal(get_cornames(names, subset = c("cor_Intercept_x", "cor_Intercept_y")), 
               c("cor(Intercept,x)", "cor(Intercept,y)"))
})

test_that("Test that get_cov_matrix returns appropriate dimensions", {
  sd <- cbind(1:10,11:20); cor <- cbind(seq(-0.5, 0.4, 0.1))
  expect_equal(dim(get_cov_matrix(sd = sd, cor = cor)$cov), c(10,2,2))
  expect_equal(dim(get_cov_matrix(sd = sd)$cov), c(10,2,2))
})

test_that("Test that evidence_ratio returns expected results", {
  post_samples <- c(-4:10); prior_samples <- c(-2:12)
  expect_true(evidence_ratio(x = post_samples, prior_samples = prior_samples) > 1)
  expect_equal(evidence_ratio(x = post_samples, cut = 0.5, wsign = "greater"), 10/5)
  expect_equal(evidence_ratio(x = post_samples, cut = 0.5, wsign = "less"), 5/10)
})

test_that("Test that expand_matrix returns expected results", {
  A <- matrix(1:6, 3, 2); x <- c(1,2,1)
  expect_equal(expand_matrix(A, x), matrix(c(1,0,3,4,0,6,0,2,0,0,5,0), 3, 4))
})

test_that("Test the find_names find all valid variable names in a string", {
  expect_equal(find_names("x + b.x - .5 + abc(a__3) : 1/2 - 0.2"), c("x", "b.x", "a__3"))
})



