  /* hurdle poisson log-PDF of a single response
   * Args:
   *   y: the response value
   *   lambda: mean parameter of the poisson distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_poisson_lpmf(int y, real lambda, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             poisson_lpmf(y | lambda) -
             log1m_exp(-lambda);
    }
  }
  /* hurdle poisson log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args:
   *   y: the response value
   *   lambda: mean parameter of the poisson distribution
   *   hu: linear predictor for hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_poisson_logit_lpmf(int y, real lambda, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             poisson_lpmf(y | lambda) -
             log1m_exp(-lambda);
    }
  }
  /* hurdle poisson log-PDF of a single response
   * log parameterization for the poisson part
   * Args:
   *   y: the response value
   *   eta: linear predictor for poisson part
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_poisson_log_lpmf(int y, real eta, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             poisson_log_lpmf(y | eta) -
             log1m_exp(-exp(eta));
    }
  }
  /* hurdle poisson log-PDF of a single response
   * log parameterization for the poisson part
   * logit parameterization of the hurdle part
   * Args:
   *   y: the response value
   *   eta: linear predictor for poisson part
   *   hu: linear predictor for hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_poisson_log_logit_lpmf(int y, real eta, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             poisson_log_lpmf(y | eta) -
             log1m_exp(-exp(eta));
    }
  }
  // hurdle poisson log-CCDF and log-CDF functions
  real hurdle_poisson_lccdf(int y, real lambda, real hu) {
    return bernoulli_lpmf(0 | hu) + poisson_lccdf(y | lambda) -
           log1m_exp(-lambda);
  }
  real hurdle_poisson_lcdf(int y, real lambda, real hu) {
    return log1m_exp(hurdle_poisson_lccdf(y | lambda, hu));
  }
