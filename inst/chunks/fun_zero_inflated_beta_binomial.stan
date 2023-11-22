  /* zero-inflated beta-binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_beta_binomial_lpmf(int y, int trials,
                                        real mu, real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         beta_binomial_lpmf(0 | trials,
                                            mu * phi,
                                            (1 - mu) * phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
    }
  }
  /* zero-inflated beta-binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_beta_binomial_logit_lpmf(int y, int trials,
                                              real mu, real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         beta_binomial_lpmf(0 | trials,
                                            mu * phi,
                                            (1 - mu) * phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
    }
  }
  // zero-inflated beta-binomial log-CCDF and log-CDF functions
  real zero_inflated_beta_binomial_lccdf(int y, int trials, real mu, real phi,
                                         real zi) {
    return bernoulli_lpmf(0 | zi) + beta_binomial_lccdf(y | trials, 
                                                        mu * phi,
                                                        (1 - mu) * phi);
  }
  real zero_inflated_beta_binomial_lcdf(int y, int trials, real mu, real phi,
                                        real zi) {
    return log1m_exp(zero_inflated_beta_binomial_lccdf(y | trials, mu, phi, zi));
  }
