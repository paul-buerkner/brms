  /* compute the cloglog link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real cloglog(real p) {
     return log(-log1m(p));
   }
