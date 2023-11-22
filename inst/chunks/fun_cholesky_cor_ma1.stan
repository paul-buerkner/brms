  /* compute the cholesky factor of a MA1 correlation matrix
   * Args:
   *   ma: MA1 autocorrelation
   *   nrows: number of rows of the covariance matrix
   * Returns:
   *   A nrows x nrows MA1 covariance matrix
   */
   matrix cholesky_cor_ma1(real ma, int nrows) {
     matrix[nrows, nrows] mat;
     mat = diag_matrix(rep_vector(1 + ma^2, nrows));
     if (nrows > 1) {
       mat[1, 2] = ma;
       for (i in 2:(nrows - 1)) {
         mat[i, i - 1] = ma;
         mat[i, i + 1] = ma;
       }
       mat[nrows, nrows - 1] = ma;
     }
     return cholesky_decompose(mat);
   }
