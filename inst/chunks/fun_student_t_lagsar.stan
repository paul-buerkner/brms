  /* student-t log-pdf for spatially lagged responses 
   * Args: 
   *   y: the response vector 
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real student_t_lagsar_lpdf(vector y, real nu, vector mu, 
                              real sigma, real rho, matrix W) {
     real K;
     matrix[rows(y), rows(y)] W_new;
     vector[rows(y)] half_pred;
     real inv_sigma2;
     K = rows(y);
     W_new = diag_matrix(rep_vector(1.0, rows(y))) - rho * W;
     half_pred  = W_new * (y - mdivide_left(W_new, mu));
     inv_sigma2 = 1 / sigma^2;
     return - K / 2 * log(nu) + lgamma((nu + K) / 2) - lgamma(nu / 2) +
            0.5 * log_determinant(crossprod(W_new) * inv_sigma2) -
            (nu + K) / 2 * log(1 + dot_self(half_pred) * inv_sigma2 / nu);
   }
