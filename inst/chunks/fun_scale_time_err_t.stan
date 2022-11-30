  /* scale and correlate time-series residuals
   * aware of time variable (for skipped time points)
   * Args:
   *   zerr: standardized and independent residuals
   *   sderr: standard deviation of the residuals
   *   chol_cor: cholesky factor of the correlation matrix
   *   n_times: number of time points in each group
   *   begin: the first observation in each group
   *   end: the last observation in each group
   *   times: vector of time points
   * Returns:
   *   vector of scaled and correlated residuals
   */
   vector scale_time_err_t(vector zerr, real sderr, matrix chol_cor,
                           int[] n_times, int[] begin, int[] end, int[] times) {
     vector[rows(zerr)] err;
     int temp_times[max(n_times)];
     for (i in 1:size(n_times)) {
       // Set up the indexing vector
       temp_times[1] = 1;
       for (j in 1:(n_times[i]-1)) {
         temp_times[j+1] = temp_times[j] +
                           (times[begin[i]+j] - times[begin[i]+j-1]);
       }
       
       err[begin[i]:end[i]] =
         sderr * chol_cor[temp_times[1:n_times[i]], temp_times[1:n_times[i]]] * 
         zerr[begin[i]:end[i]];
     }
     return err;
   }
