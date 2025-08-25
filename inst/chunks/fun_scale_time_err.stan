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
                         array[] int nobs, array[] int begin, array[] int end) {
     vector[rows(zerr)] err;
     for (i in 1:size(nobs)) {
       matrix[nobs[i], nobs[i]] L_i;
       L_i = sderr * chol_cor[1:nobs[i], 1:nobs[i]];
       err[begin[i]:end[i]] = L_i * zerr[begin[i]:end[i]];
     }
     return err;
   }
