  /* log-logistic log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: positive scale parameter
   *   shape: positive shape parameter
   * Returns:
   *   log(P(Y <= y))
   */
  real loglogistic_lcdf(real y, real mu, real shape) {
    return -log1p(pow(y / mu, -shape));
  }
  /* log-logistic log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: positive scale parameter
   *   shape: positive shape parameter
   * Returns:
   *   log(P(Y > y))
   */
  real loglogistic_lccdf(real y, real mu, real shape) {
    return -log1p(pow(y / mu, shape));
  }
