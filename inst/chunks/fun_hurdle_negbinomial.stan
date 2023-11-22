  /* hurdle negative binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_neg_binomial_lpmf(int y, real mu, real phi, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             neg_binomial_2_lpmf(y | mu, phi) -
             log1m((phi / (mu + phi))^phi);
    }
  }
  /* hurdle negative binomial log-PDF of a single response
   * logit parameterization for the hurdle part
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: phi parameter of negative binomial distribution
   *   hu: linear predictor of hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_neg_binomial_logit_lpmf(int y, real mu, real phi, real hu) {
   if (y == 0) {
     return bernoulli_logit_lpmf(1 | hu);
   } else {
     return bernoulli_logit_lpmf(0 | hu) +
            neg_binomial_2_lpmf(y | mu, phi) -
            log1m((phi / (mu + phi))^phi);
   }
  }
  /* hurdle negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi phi parameter of negative binomial distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_neg_binomial_log_lpmf(int y, real eta, real phi, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             neg_binomial_2_log_lpmf(y | eta, phi) -
             log1m((phi / (exp(eta) + phi))^phi);
    }
  }
  /* hurdle negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * logit parameterization for the hurdle part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: phi parameter of negative binomial distribution
   *   hu: linear predictor of hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_neg_binomial_log_logit_lpmf(int y, real eta, real phi, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
   } else {
      return bernoulli_logit_lpmf(0 | hu) +
             neg_binomial_2_log_lpmf(y | eta, phi) -
             log1m((phi / (exp(eta) + phi))^phi);
    }
  }
  // hurdle negative binomial log-CCDF and log-CDF functions
  real hurdle_neg_binomial_lccdf(int y, real mu, real phi, real hu) {
    return bernoulli_lpmf(0 | hu) + neg_binomial_2_lccdf(y | mu, phi) -
           log1m((phi / (mu + phi))^phi);
  }
  real hurdle_neg_binomial_lcdf(int y, real mu, real phi, real hu) {
    return log1m_exp(hurdle_neg_binomial_lccdf(y | mu, phi, hu));
  }
