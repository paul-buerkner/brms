  /* Normalizing scale factor s_k (standard deviation of unit-scale split normal)
   * Args:
   *   alpha: skewness parameter
   */
  real split_normal_s_k(real alpha) {
    real k = exp(alpha);
    return sqrt((square(k) + 1.0 / square(k) - 1.0) - (2.0 / pi()) * square(k - 1.0 / k));
  }

  /* Mode of split normal with mean mu and scale sigma
   * Args:
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   */
  real split_normal_mode(real mu, real sigma, real alpha) {
    real k = exp(alpha);
    real omega = sigma / split_normal_s_k(alpha);
    return mu - omega * (k - 1.0 / k) * sqrt(2.0 / pi());
  }

  /* split normal log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_normal_lpdf(real y, real mu, real sigma, real alpha) {
    real k = exp(alpha);
    real omega = sigma / split_normal_s_k(alpha);
    real m = split_normal_mode(mu, sigma, alpha);
    real z = (y - m) / omega;
    real u = z >= 0 ? (z / k) : (k * z);
    return log(2.0) - log(omega) - log(k + 1.0 / k) + std_normal_lpdf(u);
  }

  /* split normal log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_normal_lcdf(real y, real mu, real sigma, real alpha) {
    real k = exp(alpha);
    real omega = sigma / split_normal_s_k(alpha);
    real m = split_normal_mode(mu, sigma, alpha);
    real z = (y - m) / omega;
    real w_left = 2.0 / (1.0 + square(k));
    real w_right = (2.0 * square(k)) / (1.0 + square(k));

    if (z < 0) {
      real u = k * z;
      return log(w_left) + std_normal_lcdf(u);
    } else {
      real u = z / k;
      return log1m(w_right * exp(std_normal_lccdf(u)));
    }
  }

  /* split normal log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_normal_lccdf(real y, real mu, real sigma, real alpha) {
    real k = exp(alpha);
    real omega = sigma / split_normal_s_k(alpha);
    real m = split_normal_mode(mu, sigma, alpha);
    real z = (y - m) / omega;
    real w_left = 2.0 / (1.0 + square(k));
    real w_right = (2.0 * square(k)) / (1.0 + square(k));

    if (z >= 0) {
      real u = z / k;
      return log(w_right) + std_normal_lccdf(u);
    } else {
      real u = k * z;
      return log1m(w_left * exp(std_normal_lcdf(u)));
    }
  }
