  /* Efficient computation of the horseshoe prior
   * Args:
   *   zb: standardized population-level coefficients
   *   global: global horseshoe parameters
   *   local: local horseshoe parameters
   *   scale_global: global scale of the horseshoe prior
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector zb, vector[] local, real[] global,
                   real scale_global) {
    vector[rows(zb)] lambda;
    for (k in 1:rows(zb)) {
      lambda[k] = local[1][k] * sqrt(local[2][k]);
    }
    return zb .* lambda * global[1] * sqrt(global[2]) * scale_global;
  }
