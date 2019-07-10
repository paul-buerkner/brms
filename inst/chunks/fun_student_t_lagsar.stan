  /* student-t log-pdf for spatially lagged responses 
   * Args: 
   *   y: the response vector 
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: residual scale parameter
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   *   lambda: precomputed eigenvalues of W
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real student_t_lagsar_lpdf(vector y, real nu, vector mu, real sigma, 
                             real rho, matrix W, vector lambda) { 
    int N = rows(y);
    real K = rows(y);  // avoid integer division warning
    real inv_sigma2 = 1 / square(sigma);
    matrix[N, N] W_tilde = -rho * W;
    vector[N] half_pred;
    real log_det;
    for (n in 1:N) W_tilde[n, n] += 1;
    half_pred = W_tilde * y - mu;
    log_det = sum(log1m(rho * lambda));
    return - K / 2 * log(nu) + lgamma((nu + K) / 2) - lgamma(nu / 2) +
      0.5 * K * log(inv_sigma2) + log_det -
      (nu + K) / 2 * log(1 + dot_self(half_pred) * inv_sigma2 / nu);
  }
