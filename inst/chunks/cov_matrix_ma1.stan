  /* compute the covariance matrix for an MA1 process 
   * Args: 
   *   ma: MA1 autocorrelation 
   *   sigma: standard deviation of the MA1 process 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows MA1 covariance matrix 
   */ 
   matrix cov_matrix_ma1(real ma, real sigma, int nrows) { 
     matrix[nrows, nrows] mat; 
     mat <- diag_matrix(rep_vector(1 + ma^2, nrows)); 
     if (nrows > 1) { 
       mat[1, 2] <- ma; 
       for (i in 2:(nrows - 1)) { 
         mat[i, i - 1] <- ma; 
         mat[i, i + 1] <- ma; 
       } 
       mat[nrows, nrows - 1] <- ma; 
     } 
     return sigma^2 * mat; 
   }
