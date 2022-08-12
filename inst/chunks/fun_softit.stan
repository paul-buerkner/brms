  /* compute the softit link
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector softit(vector p) {
     return log(expm1(-p / (p - 1)));
   }
  /* compute the inverse of the sofit link
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (0, 1)
   */
   vector inv_softit(vector y) {
     return log1p_exp(y) / (1 + log1p_exp(y));
   }
