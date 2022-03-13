  /* normal log-pdf for spatially lagged residuals
   * Args:
   *   y: the response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   rho: positive autoregressive parameter
   *   W: spatial weight matrix
   *   eigenW: precomputed eigenvalues of W
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real normal_errorsar_lpdf(vector y, vector mu, real sigma,
                            real rho, data matrix W, data vector eigenW) {
    int N = rows(y);
    real inv_sigma2 = inv_square(sigma);
    matrix[N, N] W_tilde = add_diag(-rho * W, 1);
    vector[N] half_pred;
    real log_det;
    half_pred = W_tilde * (y - mu);
    log_det = sum(log1m(rho * eigenW));
    return  0.5 * N * log(inv_sigma2) + log_det -
      0.5 * dot_self(half_pred) * inv_sigma2;
  }
