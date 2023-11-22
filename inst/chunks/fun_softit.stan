  /* compute the softit link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real softit(real p) {
     return log(expm1(-p / (p - 1)));
   }
  /* compute the softit link (vectorized)
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector softit_vector(vector p) {
     return log(expm1(-p / (p - 1)));
   }
  /* compute the inverse of the sofit link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (0, 1)
   */
   real inv_softit(real y) {
     return log1p_exp(y) / (1 + log1p_exp(y));
   }
  /* compute the inverse of the sofit link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (0, 1)
   */
   vector inv_softit_vector(vector y) {
     return log1p_exp(y) / (1 + log1p_exp(y));
   }
