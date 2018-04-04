context("Tests for brmsfit helper functions")

test_that("first_greater returns expected results", {
  set.seed(12345)
  A <- cbind(1:10, 11:20, 21:30)
  x <- sample(1:25, size = 10, replace = TRUE)
  expect_equal(first_greater(A, x), c(3,3,3,3,2,1,2,2,2,3))
  expect_equal(first_greater(A, x, i = 2), c(3,3,3,3,2,2,2,2,2,3))
})

test_that("array2list performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(matrix(1:9,3,3), matrix(10:18,3,3), matrix(19:27,3,3))
  expect_equal(brms:::array2list(A), B)
})

test_that("probit and probit_approx produce similar results", {
  expect_equal(brms:::ilink(-10:10, "probit"), 
               brms:::ilink(-10:10, "probit_approx"), 
               tolerance = 1e-3)
})

test_that("ARMA covariance matrices are computed correctly", {
  ar <- 0.5
  ma <- 0.3
  sigma <- 2
  se <- sqrt(1:4)
  
  # test for AR1 cov matrix
  ar_mat <- get_cov_matrix_ar1(ar = matrix(ar), sigma = matrix(sigma), 
                               se = se, nrows = length(se))
  expected_ar_mat <- sigma^2 / (1 - ar^2) * 
                     cbind(c(1, ar, ar^2, ar^3),
                           c(ar, 1, ar, ar^2),
                           c(ar^2, ar, 1, ar),
                           c(ar^3, ar^2, ar, 1))
  expected_ar_mat <- expected_ar_mat + diag(se^2)
  expect_equal(ar_mat[1, , ], expected_ar_mat)
  
  # test for MA1 cov matrix
  ma_mat <- get_cov_matrix_ma1(ma = matrix(ma), sigma = matrix(sigma), 
                               se = se, nrows = length(se))
  expected_ma_mat <- sigma^2 * 
                     cbind(c(1+ma^2, ma, 0, 0),
                           c(ma, 1+ma^2, ma, 0),
                           c(0, ma, 1+ma^2, ma),
                           c(0, 0, ma, 1+ma^2))
  expected_ma_mat <- expected_ma_mat + diag(se^2)
  expect_equal(ma_mat[1, , ], expected_ma_mat)
  
  # test for ARMA1 cov matrix
  arma_mat <- get_cov_matrix_arma1(ar = matrix(ar), ma = matrix(ma), 
                                 sigma = matrix(sigma), 
                                 se = se, nrows = length(se))
  g0 <- 1 + ma^2 + 2 * ar * ma
  g1 <- (1 + ar * ma) * (ar + ma)
  expected_arma_mat <- sigma^2 / (1 - ar^2) * 
                       cbind(c(g0, g1, g1 * ar, g1 * ar^2),
                             c(g1, g0, g1, g1 * ar),
                             c(g1 * ar, g1, g0, g1),
                             c(g1 * ar^2, g1 * ar, g1, g0))
  expected_arma_mat <- expected_arma_mat + diag(se^2)
  expect_equal(arma_mat[1, , ], expected_arma_mat)
  
  # test for identity matrix
  ident_mat <- get_cov_matrix_ident(sigma = matrix(sigma), 
                                    se = se, nrows = length(se))
  expected_ident_mat <- diag(sigma^2 + se^2)
  expect_equal(ident_mat[1, , ], expected_ident_mat)
})

test_that("evidence_ratio returns expected results", {
  ps <- -4:10
  prs <- -2:12
  expect_true(evidence_ratio(ps, prior_samples = prs) > 1)
  expect_true(is.na(evidence_ratio(ps)))
  expect_equal(evidence_ratio(ps, cut = 0.5, wsign = "greater"), 10/5)
  expect_equal(evidence_ratio(ps, cut = 0.5, wsign = "less"), 5/10)
})

test_that("find_vars finds all valid variable names in a string", {
  string <- "x + b.x - .5 + abc(a__3) : 1/2 - 0.2"
  expect_equal(find_vars(string), c("x", "b.x", "a__3"))
})

test_that(".predictor_arma runs without errors", {
  ns <- 20
  nobs <- 30
  Y = rnorm(nobs)
  J_lag = c(1:3, 3, 3, rep(c(0:3, 3), 4), 0:3, 0)
  ar <- matrix(rnorm(ns * 3), nrow = ns, ncol = 3)
  ma <- matrix(rnorm(ns * 1), nrow = ns, ncol = 1)
  eta <- matrix(rnorm(ns * nobs), nrow = ns, ncol = nobs)
  expect_equal(.predictor_arma(eta, Y = Y, J_lag = J_lag), eta)
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ar = ar))
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ma = ma))
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ar = ar, ma = ma))
})

test_that("make_conditions works correctly", {
  conds <- make_conditions(epilepsy, c("log_Base4_c", "log_Age_c"))
  expect_equal(dim(conds), c(9, 3))
  expect_equal(conds$cond__[3], "log_Base4_c = -0.75 & log_Age_c = 0.22")
})

