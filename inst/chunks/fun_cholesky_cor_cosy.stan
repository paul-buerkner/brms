  /* compute the cholesky factor of a compound symmetry correlation matrix
   * Args:
   *   cosy: compound symmetry correlation
   *   nrows: number of rows of the covariance matrix
   * Returns:
   *   A nrows x nrows covariance matrix
   */
   matrix cholesky_cor_cosy(real cosy, int nrows) {
     matrix[nrows, nrows] mat;
     for (i in 1:nrows) {
       mat[i, i] = 1;
     }
     for (i in 2:nrows) {
       for (j in 1:(i - 1)) {
         mat[i, j] = cosy;
         mat[j, i] = mat[i, j];
       }
     }
     return cholesky_decompose(mat);
   }
