  /* scale and correlate time-series residuals
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
