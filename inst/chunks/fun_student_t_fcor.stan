  /* multi-student-t log-PDF for fixed correlation matrices
   * assuming homogoneous variances
   * Args:
   *   y: response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: scale parameter
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_fcor_hom_lpdf(vector y, real nu, vector mu, real sigma,
                               data matrix chol_cor) {
    int N = rows(chol_cor);
    matrix[N, N] Cov = multiply_lower_tri_self_transpose(sigma * chol_cor);
    return multi_student_t_lpdf(y | nu, mu, Cov);
  }
  /* multi-student-t log-PDF for fixed correlation matrices
   * assuming heterogenous variances
   * Args:
   *   y: response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: scale parameter vector
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_fcor_het_lpdf(vector y, real nu, vector mu, vector sigma,
                               data matrix chol_cor) {
    int N = rows(chol_cor);
    matrix[N, N] Cov = diag_pre_multiply(sigma, chol_cor);
    Cov = multiply_lower_tri_self_transpose(Cov);
    return multi_student_t_lpdf(y | nu, mu, Cov);
  }
