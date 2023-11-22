  /* compute the cauchit link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real cauchit(real p) {
     return tan(pi() * (p - 0.5));
   }
  /* compute the cauchit link (vectorized)
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector cauchit_vector(vector p) {
     return tan(pi() * (p - 0.5));
   }
  /* compute the inverse of the cauchit link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (0, 1)
   */
   real inv_cauchit(real y) {
     return atan(y) / pi() + 0.5;
   }
  /* compute the inverse of the cauchit link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (0, 1)
   */
   vector inv_cauchit_vector(vector y) {
     return atan(y) / pi() + 0.5;
   }
