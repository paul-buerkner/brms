  /* squareplus inverse link function (squareplus itself)
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
   *   x: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector inv_squareplus(vector x) {
     return (square(x) - 1) ./ x;
   }
