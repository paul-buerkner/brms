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
    real lp = 0.0;
    int I = size(nobs);
    matrix[rows(chol_cor), cols(chol_cor)] L = sigma * chol_cor;
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L_i = L[1:nobs[i], 1:nobs[i]];
      lp += multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L_i
      );
    }
    return lp;
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
    real lp = 0.0;
    int I = size(nobs);
    for (i in 1:I) {
      matrix[nobs[i], nobs[i]] L_i;
      L_i = diag_pre_multiply(sigma[begin[i]:end[i]], chol_cor[1:nobs[i], 1:nobs[i]]);
      lp += multi_normal_cholesky_lpdf(
        y[begin[i]:end[i]] | mu[begin[i]:end[i]], L_i
      );
    }
    return lp;
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
    real lp = 0.0;
    int I = size(nobs);
    int has_lp[I] = rep_array(0, I);
    int i = 1;
    matrix[rows(chol_cor), cols(chol_cor)] Cor;
    Cor = multiply_lower_tri_self_transpose(chol_cor);
    while (sum(has_lp) != I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      int lp_terms[I-i+1] = rep_array(0, I-i+1);
      matrix[nobs[i], nobs[i]] L_i;
      if (is_equal(iobs, sequence(1, rows(chol_cor)))) {
        // all timepoints are present in this group
        L_i = chol_cor;
      } else {
        // arbitrary subsets cannot be taken on chol_cor directly
        L_i = cholesky_decompose(Cor[iobs, iobs]);
      }
      L_i = sigma * L_i;
      has_lp[i] = 1;
      lp_terms[1] = 1;
      // find all additional groups where we have the same timepoints
      for (j in (i+1):I) {
        if (has_lp[j] == 0 && is_equal(Jtime[j], Jtime[i]) == 1) {
          has_lp[j] = 1;
          lp_terms[j-i+1] = 1;
        }
      }
      // vectorize the log likelihood by stacking the vectors
      lp += multi_normal_cholesky_lpdf(
        stack_vectors(y, nobs[i], lp_terms, begin[i:I], end[i:I]) |
        stack_vectors(mu, nobs[i], lp_terms, begin[i:I], end[i:I]), L_i
      );
      while (has_lp[i] == 1 && i != I) {
        i += 1;
      }
    }
    return lp;
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
    real lp = 0.0;
    int I = size(nobs);
    int has_lp[I] = rep_array(0, I);
    int i = 1;
    matrix[rows(chol_cor), cols(chol_cor)] Cor;
    Cor = multiply_lower_tri_self_transpose(chol_cor);
    while (sum(has_lp) != I) {
      int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
      matrix[nobs[i], nobs[i]] Lcor_i;
      matrix[nobs[i], nobs[i]] L_i;
      if (is_equal(iobs, sequence(1, rows(chol_cor)))) {
        // all timepoints are present in this group
        Lcor_i = chol_cor;
      } else {
        // arbitrary subsets cannot be taken on chol_cor directly
        Lcor_i = cholesky_decompose(Cor[iobs, iobs]);
      }
      L_i = diag_pre_multiply(sigma[begin[i]:end[i]], Lcor_i);
      lp += multi_normal_cholesky_lpdf(y[begin[i]:end[i]] | mu[begin[i]:end[i]], L_i);
      has_lp[i] = 1;
      // find all additional groups where we have the same timepoints
      for (j in (i+1):I) {
        if (has_lp[j] == 0 && is_equal(Jtime[j], Jtime[i]) == 1) {
          // group j may have different sigmas that group i
          L_i = diag_pre_multiply(sigma[begin[j]:end[j]], Lcor_i);
          lp += multi_normal_cholesky_lpdf(y[begin[j]:end[j]] | mu[begin[j]:end[j]], L_i);
          has_lp[j] = 1;
        }
      }
      while (has_lp[i] == 1 && i != I) {
        i += 1;
      }
    }
    return lp;
  }
