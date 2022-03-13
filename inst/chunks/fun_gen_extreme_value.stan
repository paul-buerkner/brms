  /* generalized extreme value log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real gen_extreme_value_lpdf(real y, real mu, real sigma, real xi) {
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return - log(sigma) - x - exp(-x);
     } else {
       real t = 1 + xi * x;
       real inv_xi = 1 / xi;
       return - log(sigma) - (1 + inv_xi) * log(t) - pow(t, -inv_xi);
     }
   }
  /* generalized extreme value log-CDF for a single response
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:
   *   log(P(Y <= y))
   */
   real gen_extreme_value_lcdf(real y, real mu, real sigma, real xi) {
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return - exp(-x);
     } else {
       return - pow(1 + xi * x, - 1 / xi);
     }
   }
  /* generalized extreme value log-CCDF for a single response
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:
   *   log(P(Y > y))
   */
   real gen_extreme_value_lccdf(real y, real mu, real sigma, real xi) {
     return log1m_exp(gen_extreme_value_lcdf(y | mu, sigma, xi));
   }
