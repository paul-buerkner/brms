  /* normal log-pdf for spatially lagged responses 
   * Args: 
   *   y: the response vector 
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real normal_lagsar_lpdf(vector y, vector mu, real sigma, 
                           real rho, matrix W) { 
    int N = rows(y);
    real inv_sigma2 = 1 / square(sigma);
    matrix[N, N] W_tilde = -rho * W;
    vector[N] half_pred;
    for (n in 1:N) W_tilde[n, n] += 1;
    half_pred = W_tilde * (y - mdivide_left(W_tilde, mu));
    return 0.5 * log_determinant(crossprod(W_tilde) * inv_sigma2) -
           0.5 * dot_self(half_pred) * inv_sigma2;
  }
