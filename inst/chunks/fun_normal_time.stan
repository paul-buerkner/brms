  /* multi-normal log-PDF for time-series covariance structures
   * assuming homogoneous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   chol_cor: cholesky factor of the correlation matrix
   *   se2: square of user defined standard errors
   *     should be set to zero if none are defined
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_hom_lpdf(vector y, vector mu, real sigma, matrix chol_cor,
                            data vector se2, int[] nobs, int[] begin, int[] end) {
    int I = size(nobs);
    int has_se = max(se2) > 0;
    vector[I] lp;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L;
      L = sigma * chol_cor[1:nobs[i], 1:nobs[i]];
      if (has_se) {
        // need to add 'se' to the covariance matrix itself
        L = multiply_lower_tri_self_transpose(L);
        L += diag_matrix(se2[begin[i]:end[i]]);
        L = cholesky_decompose(L);
      }
      lp[i] = multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
      );
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * assuming heterogenous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation vector
   *   chol_cor: cholesky factor of the correlation matrix
   *   se2: square of user defined standard errors
   *     should be set to zero if none are defined
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_het_lpdf(vector y, vector mu, vector sigma, matrix chol_cor,
                            data vector se2, int[] nobs, int[] begin, int[] end) {
    int I = size(nobs);
    int has_se = max(se2) > 0;
    vector[I] lp;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L;
      L = diag_pre_multiply(sigma[begin[i]:end[i]],
                            chol_cor[1:nobs[i], 1:nobs[i]]);
      if (has_se) {
        // need to add 'se' to the covariance matrix itself
        L = multiply_lower_tri_self_transpose(L);
        L += diag_matrix(se2[begin[i]:end[i]]);
        L = cholesky_decompose(L);
      }
      lp[i] = multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
      );
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * assuming homogoneous variances
   * allowx for flexible correlation matrix subsets
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   chol_cor: cholesky factor of the correlation matrix
   *   se2: square of user defined standard errors
   *     should be set to zero if none are defined
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_hom_flex_lpdf(vector y, vector mu, real sigma, matrix chol_cor,
                                 data vector se2, int[] nobs, int[] begin,
                                 int[] end, int[,] Jtime) {
    int I = size(nobs);
    int has_se = max(se2) > 0;
    vector[I] lp;
    matrix[rows(chol_cor), cols(chol_cor)] Cor;
    Cor = multiply_lower_tri_self_transpose(chol_cor);
    for (i in 1:I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] Cov = sigma^2 * Cor[iobs, iobs];
      if (has_se) {
        Cov += diag_matrix(se2[begin[i]:end[i]]);
      }
      lp[i] = multi_normal_lpdf(y[begin[i]:end[i]] | mu[begin[i]:end[i]], Cov);
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * assuming heterogenous variances
   * allowx for flexible correlation matrix subsets
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation vector
   *   chol_cor: cholesky factor of the correlation matrix
   *   se2: square of user defined standard errors
   *     should be set to zero if none are defined
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_het_flex_lpdf(vector y, vector mu, vector sigma, matrix chol_cor,
                                 data vector se2, int[] nobs, int[] begin,
                                 int[] end, int[,] Jtime) {
    int I = size(nobs);
    int has_se = max(se2) > 0;
    vector[I] lp;
    matrix[rows(chol_cor), cols(chol_cor)] Cor;
    Cor = multiply_lower_tri_self_transpose(chol_cor);
    for (i in 1:I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] Cov = quad_form_diag(Cor[iobs, iobs], sigma[begin[i]:end[i]]);
      if (has_se) {
        Cov += diag_matrix(se2[begin[i]:end[i]]);
      }
      lp[i] = multi_normal_lpdf(y[begin[i]:end[i]] | mu[begin[i]:end[i]], Cov);
    }
    return sum(lp);
  }
