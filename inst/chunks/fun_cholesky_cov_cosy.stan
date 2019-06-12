  /* compute the cholesky factor of a compound symmetry covariance matrix
   * Args: 
   *   cosy: compound symmetry correlation
   *   sigma: residual standard deviation of the model 
   *   nrows: number of rows of the covariance matrix 
   * Returns: 
   *   A nrows x nrows covariance matrix 
   */ 
   matrix cholesky_cov_cosy(real cosy, real sigma, int nrows) { 
     matrix[nrows, nrows] mat; 
     mat = diag_matrix(rep_vector(1, nrows)); 
     for (i in 2:nrows) { 
       for (j in 1:(i - 1)) { 
         mat[i, j] = cosy; 
         mat[j, i] = mat[i, j];
       } 
     } 
     return cholesky_decompose(sigma^2 * mat); 
   }
