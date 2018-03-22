 /* compute correlated group-level effects with 'by' variables
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: matrix of standard deviation parameters
  *   L: an array of cholesky factor correlation matrices
  *   Jby: index which grouping level belongs to which by level
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor_by(matrix z, matrix SD, matrix[] L, int[] Jby) {
    // r is stored in another dimension order than z
    matrix[cols(z), rows(z)] r;
    for (j in 1:rows(r)) {
      r[j] = (diag_pre_multiply(SD[, Jby[j]], L[Jby[j]]) * z[, j])';
    }
    return r;
  }
