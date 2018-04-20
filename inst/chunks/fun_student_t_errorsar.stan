  /* student-t log-pdf for spatially lagged residuals
   * Args: 
   *   y: the response vector 
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: residual scale parameter
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real student_t_errorsar_lpdf(vector y, real nu, vector mu, 
                               real sigma, real rho, matrix W) {
    int N = rows(y);
    real inv_sigma2 = 1 / square(sigma);
    matrix[N, N] W_tilde = -rho * W;
    vector[N] half_pred;
    for (n in 1:N) W_tilde[n, n] += 1;
    half_pred  = W_tilde * (y - mu);
    return - N / 2 * log(nu) + lgamma((nu + N) / 2) - lgamma(nu / 2) +
           0.5 * log_determinant(crossprod(W_tilde) * inv_sigma2) -
           (nu + N) / 2 * log(1 + dot_self(half_pred) * inv_sigma2 / nu);
  }
