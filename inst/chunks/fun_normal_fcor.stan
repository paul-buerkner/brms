  /* multi-normal log-PDF for fixed correlation matrices
   * assuming homogoneous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_fcor_hom_lpdf(vector y, vector mu, real sigma, data matrix chol_cor) {
    return multi_normal_cholesky_lpdf(y | mu, sigma * chol_cor);
  }
  /* multi-normal log-PDF for fixed correlation matrices
   * assuming heterogenous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation vector
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_fcor_het_lpdf(vector y, vector mu, vector sigma, data matrix chol_cor) {
    return multi_normal_cholesky_lpdf(y | mu, diag_pre_multiply(sigma, chol_cor));
  }
