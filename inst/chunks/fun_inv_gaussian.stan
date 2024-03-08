  /* inverse Gaussian log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: positive mean parameter
   *   shape: positive shape parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real inv_gaussian_lpdf(real y, real mu, real shape) {
     return 0.5 * log(shape / (2 * pi())) -
            1.5 * log(y) -
            0.5 * shape * square((y - mu) / (mu * sqrt(y)));
   }
  /* vectorized inverse Gaussian log-PDF
   * Args:
   *   y: response vector
   *   mu: positive mean parameter vector
   *   shape: positive shape parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real inv_gaussian_lpdf(vector y, vector mu, real shape) {
     return 0.5 * rows(y) * log(shape / (2 * pi())) -
            1.5 * sum(log(y)) -
            0.5 * shape * dot_self((y - mu) ./ (mu .* sqrt(y)));
   }
  /* inverse Gaussian log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: positive mean parameter
   *   shape: positive shape parameter
   * Returns:
   *   log(P(Y <= y))
   */
   real inv_gaussian_lcdf(real y, real mu, real shape) {
     return log(Phi(sqrt(shape) / sqrt(y) * (y / mu - 1)) +
              exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt(y) * (y / mu + 1)));
   }
  /* inverse Gaussian log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: positive mean parameter
   *   shape: positive shape parameter
   * Returns:
   *   log(P(Y > y))
   */
   real inv_gaussian_lccdf(real y, real mu, real shape) {
     return log1m(Phi(sqrt(shape) / sqrt(y) * (y / mu - 1)) -
              exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt(y) * (y / mu + 1)));
   }
