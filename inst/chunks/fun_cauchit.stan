  /* compute the cauchit link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real cauchit(real p) {
     return tan(pi() * (p - 0.5));
   }
  /* compute the inverse of the cauchit link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (0, 1)
   */
   real inv_cauchit(real y) {
     return cauchy_cdf(y, 0, 1);
   }
