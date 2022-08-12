  /* compute the cauchit link
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector cauchit(vector p) {
     return tan(pi() * (p - 0.5));
   }
  /* compute the inverse of the cauchit link
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (0, 1)
   */
   vector inv_cauchit(vector y) {
     return atan(y) / pi() + 0.5;
   }
