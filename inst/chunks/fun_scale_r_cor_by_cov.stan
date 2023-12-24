 /* compute correlated group-level effects with 'by' variables
  * in the presence of a within-group covariance matrix
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: matrix of standard deviation parameters
  *   L: an array of cholesky factor correlation matrices
  *   Jby: index which grouping level belongs to which by level
  *   Lcov: cholesky factor of within-group correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor_by_cov(matrix z, matrix SD, array[] matrix L,
                            array[] int Jby, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    array[size(L)] matrix[rows(L[1]), cols(L[1])] LC;
    int rows_z = rows(z);
    int rows_L = rows(L[1]);
    for (i in 1:size(LC)) {
      LC[i] = diag_pre_multiply(SD[, i], L[i]);
    }
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
              // column number of z to which z_flat[l] belongs
              int m = (l - 1) / rows_z + 1;
              r[k] = r[k] + Lcov[icov, jcov] * LC[Jby[m]][i, j] * z_flat[l];
            }
          }
        }
      }
    }
    // r is returned in another dimension order than z
    return to_matrix(r, cols(z), rows(z), 0);
  }
