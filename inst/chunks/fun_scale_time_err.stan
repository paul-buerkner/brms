  /* scale and correlate time-series residuals
   * using the Cholesky factor of the correlation matrix
   * Args:
   *   zerr: standardized and independent residuals
   *   sderr: standard deviation of the residuals
   *   chol_cor: cholesky factor of the correlation matrix
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   * Returns:
   *   vector of scaled and correlated residuals
   */
   vector scale_time_err(vector zerr, real sderr, matrix chol_cor,
                         int[] nobs, int[] begin, int[] end) {
     vector[rows(zerr)] err;
     for (i in 1:size(nobs)) {
       err[begin[i]:end[i]] =
         sderr * chol_cor[1:nobs[i], 1:nobs[i]] * zerr[begin[i]:end[i]];
     }
     return err;
   }
  /* scale and correlate time-series residuals
   * allowx for flexible correlation matrix subsets
   * Args:
   *   zerr: standardized and independent residuals
   *   sderr: standard deviation of the residuals
   *   chol_cor: cholesky factor of correlation matrix
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   *   Jtime: array of time indices per group
   * Returns:
   *   vector of scaled and correlated residuals
   */
   vector scale_time_err_flex(vector zerr, real sderr, matrix chol_cor,
                              int[] nobs, int[] begin, int[] end, int[,] Jtime) {
     vector[rows(zerr)] err;
     matrix[rows(chol_cor), cols(chol_cor)] Cor;
     Cor = multiply_lower_tri_self_transpose(chol_cor);
     for (i in 1:size(nobs)) {
       int iobs[nobs[i]] = Jtime[i, 1:nobs[i]];
       err[begin[i]:end[i]] =
         sderr * cholesky_decompose(Cor[iobs, iobs]) * zerr[begin[i]:end[i]];
     }
     return err;
   }
