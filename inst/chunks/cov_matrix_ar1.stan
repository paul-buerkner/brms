  /* compute the covariance matrix for an AR1 process 
   * Args: 
   *   ar: AR1 autocorrelation 
   *   sigma: standard deviation of the AR1 process 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows AR1 covariance matrix 
   */ 
   matrix cov_matrix_ar1(real ar, real sigma, int nrows) { 
     matrix[nrows, nrows] mat; 
     vector[nrows - 1] gamma; 
     mat <- diag_matrix(rep_vector(1, nrows)); 
     for (i in 2:nrows) { 
       gamma[i - 1] <- pow(ar, i - 1); 
       for (j in 1:(i - 1)) { 
         mat[i, j] <- gamma[i - j]; 
         mat[j, i] <- gamma[i - j]; 
       } 
     } 
     return sigma^2 / (1 - ar^2) * mat; 
   }
