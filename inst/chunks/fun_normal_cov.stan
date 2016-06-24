  /* multi-normal log-PDF for special residual covariance structures 
   * currently only ARMA effects of order 1 are implemented 
   * Args: 
   *   y: response vector 
   *   eta: linear predictor 
   *   se2: square of user defined standard errors 
   *        will be set end zero if non are defined 
   *   I: number of groups 
   *   begin: the first observation in each group 
   *   end: the last observation in each group 
   *   nobs: number of observations in each group 
   *   res_cov_matrix: AR1, MA1, or ARMA1 covariance matrix; 
   * Returns: 
   *   sum of the log-PDF values of all observations 
   */ 
   real normal_cov_lpdf(vector y, vector eta, vector se2,  
                        int I, int[] begin, int[] end, 
                        int[] nobs, matrix res_cov_matrix) { 
     vector[I] lp; 
     for (i in 1:I) { 
       matrix[nobs[i], nobs[i]] Sigma; 
       Sigma = res_cov_matrix[1:nobs[i], 1:nobs[i]] 
                + diag_matrix(se2[begin[i]:end[i]]); 
       Sigma = cholesky_decompose(Sigma); 
       lp[i] = multi_normal_cholesky_lpdf(y[begin[i]:end[i]] |
                                          eta[begin[i]:end[i]], Sigma); 
     }                        
     return sum(lp); 
   }
