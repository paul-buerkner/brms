  /* compute the logm1 link
   * Args:
   *   p: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector logm1(vector y) {
     return log(y - 1.0);
   }
  /* compute the inverse of the logm1 link
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a positive vector
   */
   vector expp1(vector y) {
     return exp(y) + 1.0;
   }
