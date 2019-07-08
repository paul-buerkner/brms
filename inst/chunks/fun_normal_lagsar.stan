  /* normal log-pdf for spatially lagged responses 
   * Args: 
   *   y: the response vector 
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   *   lambda: precomputed eigen values of W
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real normal_lagsar_lpdf(vector y, vector mu, real sigma, 
                          real rho, matrix W, vector lambda) { 
    int N = rows(y);
    real inv_sigma2 = 1 / square(sigma);
    matrix[N, N] W_tilde = -rho * W;
    vector[N] half_pred;
    real log_det;
    for (n in 1:N) W_tilde[n, n] += 1;
    half_pred = W_tilde * y - mu;
    log_det = sum(log1m(rho * lambda));
    return  0.5 * N * log(inv_sigma2) - 
      0.5 * dot_self(half_pred) * inv_sigma2 + log_det;
  }
