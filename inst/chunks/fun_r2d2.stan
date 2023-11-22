  /* compute scale parameters of the R2D2 prior
   * Args:
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   scale parameter vector of the R2D2 prior
   */
  vector scales_R2D2(vector phi, real tau2) {
    return sqrt(phi * tau2);
  }

