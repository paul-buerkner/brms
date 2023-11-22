  /* zero-inflated poisson log-PDF of a single response
   * Args:
   *   y: the response value
   *   lambda: mean parameter of the poisson distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_poisson_lpmf(int y, real lambda, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         poisson_lpmf(0 | lambda));
    } else {
      return bernoulli_lpmf(0 | zi) +
             poisson_lpmf(y | lambda);
    }
  }
  /* zero-inflated poisson log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   lambda: mean parameter of the poisson distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_poisson_logit_lpmf(int y, real lambda, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         poisson_lpmf(0 | lambda));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             poisson_lpmf(y | lambda);
    }
  }
  /* zero-inflated poisson log-PDF of a single response
   * log parameterization for the poisson part
   * Args:
   *   y: the response value
   *   eta: linear predictor for poisson distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_poisson_log_lpmf(int y, real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         poisson_log_lpmf(0 | eta));
    } else {
      return bernoulli_lpmf(0 | zi) +
             poisson_log_lpmf(y | eta);
    }
  }
  /* zero-inflated poisson log-PDF of a single response
   * log parameterization for the poisson part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   eta: linear predictor for poisson distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_poisson_log_logit_lpmf(int y, real eta, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         poisson_log_lpmf(0 | eta));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             poisson_log_lpmf(y | eta);
    }
  }
  // zero-inflated poisson log-CCDF and log-CDF functions
  real zero_inflated_poisson_lccdf(int y, real lambda, real zi) {
    return bernoulli_lpmf(0 | zi) + poisson_lccdf(y | lambda);
  }
  real zero_inflated_poisson_lcdf(int y, real lambda, real zi) {
    return log1m_exp(zero_inflated_poisson_lccdf(y | lambda, zi));
  }
