  /* multi-student-t log-PDF for special residual covariance structures 
   * Args: 
   *   y: response vector 
   *   nu: degrees of freedom parameter 
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
   real student_t_cov_lpdf(vector y, real nu, vector mu, matrix chol_cov, 
                           vector se2, int[] nobs, int[] begin, int[] end) { 
     int I = size(nobs);
     vector[I] lp; 
     for (i in 1:I) { 
       matrix[nobs[i], nobs[i]] Cov = chol_cov[1:nobs[i], 1:nobs[i]];
       Cov = multiply_lower_tri_self_transpose(Cov) 
         + diag_matrix(se2[begin[i]:end[i]]);
       lp[i] = multi_student_t_lpdf(
          y[begin[i]:end[i]] | nu, mu[begin[i]:end[i]], Cov
       ); 
     }                        
     return sum(lp); 
   }
