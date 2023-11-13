  /* multi-student-t log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming homogoneous variances
   * and known standard errors
   * Args:
   *   y: response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: scale parameter
   *   se2: square of user defined standard errors
   *   chol_cor: cholesky factor of the correlation matrix
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_time_hom_se_lpdf(vector y, real nu, vector mu, real sigma,
                                  data vector se2, matrix chol_cor, array[] int nobs,
                                  array[] int begin, array[] int end) {
    int I = size(nobs);
    vector[I] lp;
    matrix[rows(chol_cor), cols(chol_cor)] Cov;
    Cov = multiply_lower_tri_self_transpose(sigma * chol_cor);
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] Cov_i = Cov[1:nobs[i], 1:nobs[i]];
      // need to add 'se' to the covariance matrix itself
      Cov_i += diag_matrix(se2[begin[i]:end[i]]);
      lp[i] = multi_student_t_lpdf(y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov_i);
    }
    return sum(lp);
  }
  /* multi-student-t log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming heterogenous variances
   * and known standard errors
   * Deviating Args:
   *   sigma: scale parameter vector
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_time_het_se_lpdf(vector y, real nu, vector mu, vector sigma,
                                  data vector se2, matrix chol_cor, array[] int nobs,
                                  array[] int begin, array[] int end) {
    int I = size(nobs);
    vector[I] lp;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] Cov_i;
      Cov_i = diag_pre_multiply(sigma[begin[i]:end[i]], chol_cor[1:nobs[i], 1:nobs[i]]);
      Cov_i = multiply_lower_tri_self_transpose(Cov_i);
      Cov_i += diag_matrix(se2[begin[i]:end[i]]);
      lp[i] = multi_student_t_lpdf(y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov_i);
    }
    return sum(lp);
  }
  /* multi-student-t log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming homogoneous variances
   * and known standard errors
   * allows for flexible correlation matrix subsets
   * Deviating Args:
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_time_hom_se_flex_lpdf(vector y, real nu, vector mu, real sigma,
                                       data vector se2, matrix chol_cor, array[] int nobs,
                                       array[] int begin, array[] int end, array[,] int Jtime) {
    int I = size(nobs);
    vector[I] lp;
    matrix[rows(chol_cor), cols(chol_cor)] Cov;
    Cov = multiply_lower_tri_self_transpose(sigma * chol_cor);
    for (i in 1:I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] Cov_i = Cov[iobs, iobs];
      Cov_i += diag_matrix(se2[begin[i]:end[i]]);
      lp[i] = multi_student_t_lpdf(y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov_i);
    }
    return sum(lp);
  }
  /* multi-student-t log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming heterogenous variances
   * and known standard errors
   * allows for flexible correlation matrix subsets
   * Deviating Args:
   *   sigma: scale parameter vector
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_time_het_se_flex_lpdf(vector y, real nu, vector mu, vector sigma,
                                       data vector se2, matrix chol_cor, array[] int nobs,
                                       array[] int begin, array[] int end, array[,] int Jtime) {
    int I = size(nobs);
    vector[I] lp;
    matrix[rows(chol_cor), cols(chol_cor)] Cor;
    Cor = multiply_lower_tri_self_transpose(chol_cor);
    for (i in 1:I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] Cov_i;
      Cov_i = quad_form_diag(Cor[iobs, iobs], sigma[begin[i]:end[i]]);
      Cov_i += diag_matrix(se2[begin[i]:end[i]]);
      lp[i] = multi_student_t_lpdf(y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov_i);
    }
    return sum(lp);
  }
