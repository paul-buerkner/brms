  /* compute the softit link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real softit(real p) {
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
