  /* zero-inflated negative binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_lpmf(int y, real mu, real phi,
                                       real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_logit_lpmf(int y, real mu,
                                             real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_lpmf(int y, real eta,
                                           real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_logit_lpmf(int y, real eta,
                                                 real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  // zero_inflated negative binomial log-CCDF and log-CDF functions
  real zero_inflated_neg_binomial_lccdf(int y, real mu, real phi, real hu) {
    return bernoulli_lpmf(0 | hu) + neg_binomial_2_lccdf(y | mu, phi);
  }
  real zero_inflated_neg_binomial_lcdf(int y, real mu, real phi, real hu) {
    return log1m_exp(zero_inflated_neg_binomial_lccdf(y | mu, phi, hu));
  }
