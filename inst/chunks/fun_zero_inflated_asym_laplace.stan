  /* zero-inflated asymmetric laplace log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_asym_laplace_lpdf(real y, real mu, real sigma,
                                       real quantile, real zi) {
    if (y == 0) {
      return bernoulli_lpmf(1 | zi);
    } else {
      return bernoulli_lpmf(0 | zi) +
             asym_laplace_lpdf(y | mu, sigma, quantile);
    }
  }
  /* zero-inflated asymmetric laplace log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   *   zi: linear predictor of the zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_asym_laplace_logit_lpdf(real y, real mu, real sigma,
                                             real quantile, real zi) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | zi);
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             asym_laplace_lpdf(y | mu, sigma, quantile);
    }
  }
  // zero-inflated asymmetric laplace log-CDF function
  real zero_inflated_asym_laplace_lcdf(real y, real mu, real sigma,
                                       real quantile, real zi) {
    if (y < 0) {
      return bernoulli_lpmf(0 | zi) +
             asym_laplace_lcdf(y | mu, sigma, quantile);
    } else {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         asym_laplace_lcdf(y | mu, sigma, quantile));
    }
  }
  // zero-inflated asymmetric laplace log-CCDF function
  real zero_inflated_asym_laplace_lccdf(real y, real mu, real sigma,
                                        real quantile, real zi) {
    if (y > 0) {
      return bernoulli_lpmf(0 | zi) +
             asym_laplace_lccdf(y | mu, sigma, quantile);
    } else {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         asym_laplace_lccdf(y | mu, sigma, quantile));
    }
  }
