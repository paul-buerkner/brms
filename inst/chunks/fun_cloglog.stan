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
   vector cloglog(vector p) {
     return log(-log1m(p));
   }
  /* compute the log of the inverse of the cloglog link, which is y to
   * double precision below y = -37 and 0 from y = 709, just below where
   * exp(y) overflows and would make gradients NaN
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-Inf, 0)
   */
   real log_inv_cloglog(real y) {
     if (y < -37) {
       return y;
     } else if (y < 709) {
       return log1m_exp(-exp(y));
     }
     return 0;
   }
  /* compute the log of one minus the inverse of the cloglog link, which is
   * -Inf where exp(y) overflows rather than -exp(y), whose infinite gradient
   * would make others NaN
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-Inf, 0)
   */
   real log1m_inv_cloglog(real y) {
     if (y > 709.782712893384) {
       return negative_infinity();
     }
     return -exp(y);
   }
