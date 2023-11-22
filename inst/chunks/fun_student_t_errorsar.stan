  /* student-t log-pdf for spatially lagged residuals
   * Args:
   *   y: the response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: residual scale parameter
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   *   eigenW: precomputed eigenvalues of W
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real student_t_errorsar_lpdf(vector y, real nu, vector mu, real sigma,
                               real rho, data matrix W, data vector eigenW) {
    int N = rows(y);
    real K = rows(y);  // avoid integer division warning
    real inv_sigma2 = inv_square(sigma);
    matrix[N, N] W_tilde = add_diag(-rho * W, 1);
    vector[N] half_pred;
    real log_det;
    half_pred = W_tilde * (y - mu);
    log_det = sum(log1m(rho * eigenW));
    return - K / 2 * log(nu) + lgamma((nu + K) / 2) - lgamma(nu / 2) +
      0.5 * K * log(inv_sigma2) + log_det -
      (nu + K) / 2 * log1p(dot_self(half_pred) * inv_sigma2 / nu);
  }
