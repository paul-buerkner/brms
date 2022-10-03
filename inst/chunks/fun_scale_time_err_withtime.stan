  /* scale and correlate time-series residuals
   * aware of time variable (for skipped time points)
   * Args:
   *   zerr: standardized and independent residuals
   *   sderr: standard deviation of the residuals
   *   chol_cor: cholesky factor of the correlation matrix
   *   nobs: number of observations in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   *   times: 
   * Returns:
   *   vector of scaled and correlated residuals
   */
   vector scale_time_err_t(vector zerr, real sderr, matrix chol_cor,
                           int[] nobs, int[] begin, int[] end, int[] times) {
     vector[rows(zerr)] err;
     int temp_times[max(begin - end)];
     for (i in 1:size(nobs)) {
       // Set up the indexing vector
       temp_times[1] = 1;
       for (j in 1:(end[i]-begin[i]-1)) {
         temp_times[j+1] = temp_times[j] +
                           (times[begin[i] + j + 1] - 
                           times[begin[i] + j])
       }
       
       err[begin[i]:end[i]] =
         sderr * chol_cor[temp_times[1:nobs[i]], temp_times[1:nobs[i]]] * 
         zerr[begin[i]:end[i]];
     }
     return err;
   }
