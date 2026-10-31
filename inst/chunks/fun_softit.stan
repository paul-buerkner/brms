  /* compute the softit link
   * Args:
   *   p: a scalar in (0, 1)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real softit(real p) {
     return log(expm1(-p / (p - 1)));
   }
  /* compute the softit link (vectorized)
   * Args:
   *   p: a vector in (0, 1)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector softit(vector p) {
     return log(expm1(-p ./ (p - 1)));
   }
  /* compute the inverse of the sofit link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (0, 1)
   */
   real inv_softit(real y) {
     return log1p_exp(y) / (1 + log1p_exp(y));
   }
  /* compute the inverse of the sofit link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (0, 1)
   */
   vector inv_softit(vector y) {
     return log1p_exp(y) ./ (1 + log1p_exp(y));
   }
  /* compute the log of the inverse of the softit link, which is y to
   * double precision below y = -37
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-Inf, 0)
   */
   real log_inv_softit(real y) {
     if (y < -37) {
       return y;
     }
     return -log1p(inv(log1p_exp(y)));
   }
  /* compute the log of one minus the inverse of the softit link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-Inf, 0)
   */
   real log1m_inv_softit(real y) {
     return -log1p(log1p_exp(y));
   }
