  /* compute the cholesky factor of an ARMA1 covariance matrix
   * Args: 
   *   ar: AR1 autocorrelation 
   *   ma: MA1 autocorrelation 
   *   sigma: standard deviation of the ARMA1 process 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows matrix 
   */ 
   matrix cholesky_cov_arma1(real ar, real ma, real sigma, int nrows) { 
     matrix[nrows, nrows] mat; 
     vector[nrows] gamma; 
     mat = diag_matrix(rep_vector(1 + ma^2 + 2 * ar * ma, nrows)); 
     gamma[1] = (1 + ar * ma) * (ar + ma); 
     for (i in 2:nrows) { 
       gamma[i] = gamma[1] * pow(ar, i - 1); 
       for (j in 1:(i - 1)) { 
         mat[i, j] = gamma[i - j]; 
         mat[j, i] = gamma[i - j]; 
       } 
     } 
     return cholesky_decompose(sigma^2 / (1 - ar^2) * mat); 
   }
