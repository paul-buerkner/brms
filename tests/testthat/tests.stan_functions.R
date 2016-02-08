test_that("self-defined Stan functions work correctly", {
  rstan::expose_stan_functions(brms:::new_stan_functions)
  # arma covariance matrices 
  cov_ar1_R <- get_cov_matrix_ar1(ar = matrix(0.5), sigma = 2, 
                                  sq_se = 0, nrows = 3)[1, , ]
  expect_equal(cov_matrix_ar1(0.5, 2, 3), cov_ar1_R)
  cov_ma1_R <- matrix(get_cov_matrix_ma1(ma = matrix(-0.3), sigma = 3, 
                                         sq_se = 0, nrows = 1)[1, , ])
  expect_equal(cov_matrix_ma1(-0.3, 3, 1), cov_ma1_R)
  cov_arma1_R <- get_cov_matrix_arma1(ar = matrix(-0.5), ma = matrix(0.7), 
                                      sigma = 4, sq_se = 0, nrows = 5)[1, , ]
  expect_equal(cov_matrix_arma1(-0.5, 0.7, 4, 5), cov_arma1_R)
  # kronecker product
  A <- matrix(c(3, 2, 1, 2, 4, 1, 1, 1, 5), nrow = 3)
  B <- matrix(c(3, 2, 2, 4), nrow = 2)
  sd <- c(2, 7)
  expect_equal(t(kronecker_cholesky(A, t(chol(B)), sd)),
               chol(kronecker(A, diag(sd) %*% B %*% diag(sd))))
  # cauchit link
  expect_equal(inv_cauchit(1.5), pcauchy(1.5))
  # 
})
