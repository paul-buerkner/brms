  /* squareplus inverse link function (squareplus itself)
   * Args:
   *   x: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
   real squareplus(real x) {
     return (x + sqrt(square(x) + 4)) / 2;
   }
  /* squareplus inverse link function (squareplus itself; vectorized)
   * Args:
   *   x: a vector in (-Inf, Inf)
   * Returns:
   *   a positive vector
   */
   vector squareplus(vector x) {
     return (x + sqrt(square(x) + 4)) / 2;
   }
  /* squareplus link function (inverse squareplus)
   * Args:
   *   x: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real inv_squareplus(real x) {
     return (square(x) - 1) ./ x;
   }
  /* squareplus link function (inverse squareplus; vectorized)
   * Args:
   *   x: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector inv_squareplus(vector x) {
     return (square(x) - 1) ./ x;
   }
