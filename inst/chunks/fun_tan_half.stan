  /* compute the tan_half link
   * Args:
   *   x: a scalar in (-pi, pi)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real tan_half(real x) {
     return tan(x / 2);
   }
  /* compute the inverse of the tan_half link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-pi, pi)
   */
   real inv_tan_half(real y) {
     return 2 * atan(y);
   }
