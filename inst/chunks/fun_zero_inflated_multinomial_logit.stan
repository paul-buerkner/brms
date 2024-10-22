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
  /* Zero-inflated-dirichlet-multinomial-logit log-PDF
   * Args:
   *   y: vector of real response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real zero_inflated_multinomial_logit2_lpmf(array[] int y, vector mu, vector zi, real phi) {
      real ptarget = 0.0;
      vector [size(mu)] lambda = rep(0.0,size(mu));
      alpha = softmax(mu) ;
      int sum_y = sum(y);
      for (i in 1:size(mu)) {
        ptarget += zero_inflated_poisson_logit_lpmf(y[i]| alpha[i] * sum_y,zi[i]);
      }
      return ptarget;
  }
