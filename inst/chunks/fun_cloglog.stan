  /* compute the cloglog link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real cloglog(real p) {
     return log(-log1m(p));
   }
  /* compute the cloglog link (vectorized)
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector cloglog_vector(vector p) {
     return log(-log1m(p));
   }
