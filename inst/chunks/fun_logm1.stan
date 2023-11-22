  /* compute the logm1 link
   * Args:
   *   p: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real logm1(real y) {
     return log(y - 1.0);
   }
  /* compute the logm1 link (vectorized)
   * Args:
   *   p: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector logm1_vector(vector y) {
     return log(y - 1.0);
   }
  /* compute the inverse of the logm1 link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
   real expp1(real y) {
     return exp(y) + 1.0;
   }
  /* compute the inverse of the logm1 link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a positive vector
   */
   vector expp1_vector(vector y) {
     return exp(y) + 1.0;
   }
