  /* compute the cloglog link
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector cloglog(vector p) {
     return log(-log1m(p));
   }
