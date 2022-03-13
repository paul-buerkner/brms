  /* zero-inflated binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   theta: probability parameter of the binomial part
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_lpmf(int y, int trials,
                                   real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   theta: probability parameter of the binomial part
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_logit_lpmf(int y, int trials,
                                         real theta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_lpmf(0 | trials, theta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_lpmf(y | trials, theta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the binomial part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   eta: linear predictor for binomial part
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_blogit_lpmf(int y, int trials,
                                          real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  /* zero-inflated binomial log-PDF of a single response
   * logit parameterization of the binomial part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   trials: number of trials of the binomial part
   *   eta: linear predictor for binomial part
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_binomial_blogit_logit_lpmf(int y, int trials,
                                                real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         binomial_logit_lpmf(0 | trials, eta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             binomial_logit_lpmf(y | trials, eta);
    }
  }
  // zero-inflated binomial log-CCDF and log-CDF functions
  real zero_inflated_binomial_lccdf(int y, int trials, real theta, real zi) {
    return bernoulli_lpmf(0 | zi) + binomial_lccdf(y | trials, theta);
  }
  real zero_inflated_binomial_lcdf(int y, int trials, real theta, real zi) {
    return log1m_exp(zero_inflated_binomial_lccdf(y | trials, theta, zi));
  }
