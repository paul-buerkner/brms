  /* multi-normal log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming homogoneous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   chol_cor: cholesky factor of the correlation matrix
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_hom_lpdf(vector y, vector mu, real sigma, matrix chol_cor,
                            int[] nobs, int[] begin, int[] end) {
    int I = size(nobs);
    vector[I] lp;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L;
      L = sigma * chol_cor[1:nobs[i], 1:nobs[i]];
      lp[i] = multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
      );
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming heterogenous variances
   * Deviating Args:
   *   sigma: residual standard deviation vector
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_het_lpdf(vector y, vector mu, vector sigma, matrix chol_cor,
                            int[] nobs, int[] begin, int[] end) {
    int I = size(nobs);
    vector[I] lp;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L;
      L = diag_pre_multiply(sigma[begin[i]:end[i]], chol_cor[1:nobs[i], 1:nobs[i]]);
      lp[i] = multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
      );
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming homogoneous variances
   * allows for flexible correlation matrix subsets
   * Deviating Args:
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_hom_flex_lpdf(vector y, vector mu, real sigma, matrix chol_cor,
                                 int[] nobs, int[] begin, int[] end, int[,] Jtime) {
    int I = size(nobs);
    vector[I] lp;
    int has_lp[I] = rep_array(0, I);
    int i = 1;
    while (sum(has_lp) != I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] L = diag_pre_multiply(rep_vector(sigma, nobs[i]), chol_cor[iobs, iobs]);
      lp[i] = multi_normal_cholesky_lpdf(y[begin[i]:end[i]] | mu[begin[i]:end[i]], L);
      has_lp[i] = 1;
      // find all additional groups where we have the same timepoints
      for (j in (i+1):I) {
        if (has_lp[j] == 0 && is_equal(Jtime[j], Jtime[i]) == 1) {
          lp[j] = multi_normal_cholesky_lpdf(y[begin[j]:end[j]] | mu[begin[j]:end[j]], L);
          has_lp[j] = 1;
        }
      }
      while (has_lp[i] == 1 && i != I) {
        i += 1;
      }
    }
    return sum(lp);
  }
  /* multi-normal log-PDF for time-series covariance structures
   * in Cholesky parameterization and assuming heterogenous variances
   * allows for flexible correlation matrix subsets
   * Deviating Args:
   *   sigma: residual standard deviation vectors
   *   Jtime: array of time indices per group
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_time_het_flex_lpdf(vector y, vector mu, vector sigma, matrix chol_cor,
                                 int[] nobs, int[] begin, int[] end, int[,] Jtime) {
    int I = size(nobs);
    vector[I] lp;
    int has_lp[I] = rep_array(0, I);
    int i = 1;
    while (sum(has_lp) != I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] L = diag_pre_multiply(sigma[begin[i]:end[i]], chol_cor[iobs, iobs]);
      lp[i] = multi_normal_cholesky_lpdf(y[begin[i]:end[i]] | mu[begin[i]:end[i]], L);
      has_lp[i] = 1;
      // find all additional groups where we have the same timepoints
      for (j in (i+1):I) {
        if (has_lp[j] == 0 && is_equal(Jtime[j], Jtime[i]) == 1) {
          lp[j] = multi_normal_cholesky_lpdf(y[begin[j]:end[j]] | mu[begin[j]:end[j]], L);
          has_lp[j] = 1;
        }
      }
      while (has_lp[i] == 1 && i != I) {
        i += 1;
      }
    }
    return sum(lp);
  }
