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

test_that("Test that ARMA covariance matrices are computed correctly", {
  ar <- 0.5
  ma <- 0.3
  sigma <- 2
  sq_se <- 1:4
  # test for AR1 cov matrix
  ar_mat <- get_cov_matrix_ar1(ar = matrix(ar), sigma = matrix(sigma), 
                               sq_se = sq_se, nrows = length(sq_se))
  expected_ar_mat <- sigma^2 / (1 - ar^2) * 
                     cbind(c(1, ar, ar^2, ar^3),
                           c(ar, 1, ar, ar^2),
                           c(ar^2, ar, 1, ar),
                           c(ar^3, ar^2, ar, 1))
  expected_ar_mat <- expected_ar_mat + diag(sq_se)
  expect_equal(ar_mat[1, , ], expected_ar_mat)
  # test for MA1 cov matrix
  ma_mat <- get_cov_matrix_ma1(ma = matrix(ma), sigma = matrix(sigma), 
                               sq_se = sq_se, nrows = length(sq_se))
  expected_ma_mat <- sigma^2 * 
                     cbind(c(1+ma^2, ma, 0, 0),
                           c(ma, 1+ma^2, ma, 0),
                           c(0, ma, 1+ma^2, ma),
                           c(0, 0, ma, 1+ma^2))
  expected_ma_mat <- expected_ma_mat + diag(sq_se)
  expect_equal(ma_mat[1, , ], expected_ma_mat)
  # test for ARMA1 cov matrix
  arma_mat <- get_cov_matrix_arma1(ar = matrix(ar), ma = matrix(ma), 
                                 sigma = matrix(sigma), 
                                 sq_se = sq_se, nrows = length(sq_se))
  g0 <- 1 + ma^2 + 2 * ar * ma
  g1 <- (1 + ar * ma) * (ar + ma)
  expected_arma_mat <- sigma^2 / (1 - ar^2) * 
                       cbind(c(g0, g1, g1 * ar, g1 * ar^2),
                             c(g1, g0, g1, g1 * ar),
                             c(g1 * ar, g1, g0, g1),
                             c(g1 * ar^2, g1 * ar, g1, g0))
  expected_arma_mat <- expected_arma_mat + diag(sq_se)
  expect_equal(arma_mat[1, , ], expected_arma_mat)
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

test_that("Test that get_summary returns corrects dims and names", {
  col_names <- c("Estimate", "Est.Error", "5%ile", "50%ile", "95%ile")
  samples_2dim <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  summary_2dim <- get_summary(samples_2dim, probs = c(0.05, 0.5, 0.95))
  expect_equal(dim(summary_2dim), c(10, 5))
  expect_equal(dimnames(summary_2dim), 
               list(as.character(1:10), col_names))
  
  col_names <- c("Estimate", "Est.Error", "2.5%ile", "97.5%ile")
  samples_3dim <- array(sample(1:10, x = 4000, replace = TRUE), 
                        dim = c(100, 10, 4))
  summary_3dim <- get_summary(samples_3dim)
  expect_equal(dim(summary_3dim), c(10, 4, 4))
  expect_equal(dimnames(summary_3dim), 
               list(as.character(1:10), col_names, 
                    paste0("P(Y = ", 1:4, ")")))
  expect_error(get_summary(rnorm(100)))
})

