  /* scale and correlate residuals 
   * Args: 
   *   zerr: standardized and independent residuals
   *   chol_cov: cholesky factor of the covariance matrix
   *   nobs: number of observations in each group 
   *   begin: the first observation in each group 
   *   end: the last observation in each group 
   * Returns: 
   *   vector of scaled and correlated residuals
   */ 
   vector scale_cov_err(vector zerr, matrix chol_cov, int[] nobs,
                        int[] begin, int[] end) { 
     vector[rows(zerr)] err; 
     for (i in 1:size(nobs)) { 
       err[begin[i]:end[i]] = 
         chol_cov[1:nobs[i], 1:nobs[i]] * zerr[begin[i]:end[i]];
     }                        
     return err; 
   }
