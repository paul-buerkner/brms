 /* compute correlated group-level effects
  * in the presence of a within-group covariance matrix
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  *   Lcov: cholesky factor of within-group correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor_cov(matrix z, vector SD, matrix L, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    matrix[rows(L), cols(L)] LC = diag_pre_multiply(SD, L);
    int rows_z = rows(z);
    int rows_L = rows(L);
    // kronecker product of cholesky factors times a vector
    for (icov in 1:rows(Lcov)) {
      for (jcov in 1:icov) {
        if (Lcov[icov, jcov] > 1e-10) {
          // avoid calculating products between unrelated individuals
          for (i in 1:rows_L) {
            for (j in 1:i) {
              // incremented element of the output vector
              int k = (rows_L * (icov - 1)) + i;
              // applied element of the input vector
              int l = (rows_L * (jcov - 1)) + j;
              r[k] = r[k] + Lcov[icov, jcov] * LC[i, j] * z_flat[l];
            }
          }
        }
      }
    }
    // r is returned in another dimension order than z
    return to_matrix(r, cols(z), rows(z), 0);
  }
