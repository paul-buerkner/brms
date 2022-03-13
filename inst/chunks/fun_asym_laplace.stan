  /* helper function for asym_laplace_lpdf
   * Args:
   *   y: the response value
   *   quantile: quantile parameter in (0, 1)
   */
   real rho_quantile(real y, real quantile) {
     if (y < 0) {
       return y * (quantile - 1);
     } else {
       return y * quantile;
     }
   }
  /* asymmetric laplace log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lpdf(real y, real mu, real sigma, real quantile) {
     return log(quantile * (1 - quantile)) -
            log(sigma) -
            rho_quantile((y - mu) / sigma, quantile);
   }
  /* asymmetric laplace log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lcdf(real y, real mu, real sigma, real quantile) {
     if (y < mu) {
       return log(quantile) + (1 - quantile) * (y - mu) / sigma;
     } else {
       return log1m((1 - quantile) * exp(-quantile * (y - mu) / sigma));
     }
   }
  /* asymmetric laplace log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lccdf(real y, real mu, real sigma, real quantile) {
     if (y < mu) {
       return log1m(quantile * exp((1 - quantile) * (y - mu) / sigma));
     } else {
       return log1m(quantile) - quantile * (y - mu) / sigma;
     }
   }
