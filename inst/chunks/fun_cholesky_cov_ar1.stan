  /* compute the cholesky factor of an AR1 covariance matrix
   * Args: 
   *   ar: AR1 autocorrelation 
   *   sigma: standard deviation of the AR1 process 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows matrix 
   */ 
   matrix cholesky_cov_ar1(real ar, real sigma, int nrows) { 
     matrix[nrows, nrows] mat; 
     vector[nrows - 1] gamma; 
     mat = diag_matrix(rep_vector(1, nrows)); 
     for (i in 2:nrows) { 
       gamma[i - 1] = pow(ar, i - 1); 
       for (j in 1:(i - 1)) { 
         mat[i, j] = gamma[i - j]; 
         mat[j, i] = gamma[i - j]; 
       } 
     } 
     return cholesky_decompose(sigma^2 / (1 - ar^2) * mat); 
   }
