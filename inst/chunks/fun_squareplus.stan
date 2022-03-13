  /* squareplus inverse link function (squareplus itself)
   * Args:
   *   x: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
   real squareplus(real x) {
     return (x + sqrt(x^2 + 4)) / 2;
   }
  /* squareplus link function (inverse squareplus)
   * Args:
   *   x: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real inv_squareplus(real x) {
     return (x^2 - 1) / x;
   }
