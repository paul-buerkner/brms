  /* zero-inflated beta log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real zero_inflated_beta_lpdf(real y, real mu, real phi, real zi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) {
       return bernoulli_lpmf(1 | zi);
     } else {
       return bernoulli_lpmf(0 | zi) +
              beta_lpdf(y | shape[1], shape[2]);
     }
   }
  /* zero-inflated beta log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real zero_inflated_beta_logit_lpdf(real y, real mu, real phi, real zi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) {
       return bernoulli_logit_lpmf(1 | zi);
     } else {
       return bernoulli_logit_lpmf(0 | zi) +
              beta_lpdf(y | shape[1], shape[2]);
     }
   }
  // zero-inflated beta log-CCDF and log-CDF functions
  real zero_inflated_beta_lccdf(real y, real mu, real phi, real zi) {
    row_vector[2] shape = [mu * phi, (1 - mu) * phi];
    return bernoulli_lpmf(0 | zi) + beta_lccdf(y | shape[1], shape[2]);
  }
  real zero_inflated_beta_lcdf(real y, real mu, real phi, real zi) {
    return log1m_exp(zero_inflated_beta_lccdf(y | mu, phi, zi));
  }
