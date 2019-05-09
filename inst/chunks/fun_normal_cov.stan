  /* multi-normal log-PDF for special residual covariance structures 
   * Args: 
   *   y: response vector 
   *   mu: mean parameter vector
   *   chol_cov: cholesky factor of the covariance matrix
   *   se2: square of user defined standard errors 
   *     should be set to zero if none are defined 
   *   nobs: number of observations in each group 
   *   begin: the first observation in each group 
   *   end: the last observation in each group 
   * Returns: 
   *   sum of the log-PDF values of all observations 
   */ 
   real normal_cov_lpdf(vector y, vector mu, matrix chol_cov, vector se2,  
                        int[] nobs, int[] begin, int[] end) {
     int I = size(nobs);
     int has_se = max(se2) > 0;
     vector[I] lp; 
     for (i in 1:I) { 
       matrix[nobs[i], nobs[i]] L = chol_cov[1:nobs[i], 1:nobs[i]];
       if (has_se) {
         // need to add 'se' to the correlation matrix itself
         L = multiply_lower_tri_self_transpose(L) 
           + diag_matrix(se2[begin[i]:end[i]]);
         L = cholesky_decompose(L);
       }
       lp[i] = multi_normal_cholesky_lpdf(
         y[begin[i]:end[i]] | mu[begin[i]:end[i]], L
       ); 
     }                        
     return sum(lp); 
   }
