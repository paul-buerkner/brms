  /* First absolute moment of Student-t with nu degrees of freedom: E[|T_nu|]
   * Args:
   *   nu: degrees of freedom parameter (> 1)
   */
  real fs_abs_mean(real nu) {
    return 2 * sqrt(nu) * exp(lgamma((nu + 1) / 2) - lgamma(nu / 2)) / ((nu - 1) * sqrt(pi()));
  }

  /* Normalizing scale factor s_k (asymptotic SD of unit-scale split normal)
   * Args:
   *   alpha: skewness parameter
   */
  real split_s_k(real alpha) {
    real k = exp(alpha);
    return sqrt((square(k) + 1 / square(k) - 1) - (2 / pi()) * square(k - 1 / k));
  }

  /* Mode of split Student-t with mean mu and scale sigma
   * Args:
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   *   nu: degrees of freedom parameter (> 1)
   */
  real split_mode(real mu, real sigma, real alpha, real nu) {
    real k = exp(alpha);
    real omega = sigma / split_s_k(alpha);
    return mu - omega * (k - 1 / k) * fs_abs_mean(nu);
  }

  /* split Student-t log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   *   nu: degrees of freedom parameter (> 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_student_lpdf(real y, real mu, real sigma, real alpha, real nu) {
    real k = exp(alpha);
    real omega = sigma / split_s_k(alpha);
    real m = split_mode(mu, sigma, alpha, nu);
    real z = (y - m) / omega;
    real u = z >= 0 ? (z / k) : (k * z);
    return log(2.0) - log(omega) - log(k + 1 / k) + student_t_lpdf(u | nu, 0, 1);
  }

  /* split Student-t log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   *   nu: degrees of freedom parameter (> 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_student_lcdf(real y, real mu, real sigma, real alpha, real nu) {
    real k = exp(alpha); 
    real omega = sigma / split_s_k(alpha);
    real m = split_mode(mu, sigma, alpha, nu);
    real z = (y - m) / omega;
    real w_left  = 2.0 / (1.0 + square(k));
    real w_right = 2.0 * square(k) / (1.0 + square(k));
    
    if (z < 0) {
      real u = k * z;
      // Guard at u = 0 ensures finite gradient for Stan's autodiff at the mode
      real log_T = (u == 0) ? log(0.5) + 2 * u * exp(student_t_lpdf(0 | nu, 0, 1))
                            : student_t_lcdf(u | nu, 0, 1);
      return log(w_left) + log_T;
    } else {
      real u = z / k;
      real log_1m_T = (u == 0) ? log(0.5) - 2 * u * exp(student_t_lpdf(0 | nu, 0, 1))
                                : student_t_lccdf(u | nu, 0, 1);
      return log1m(w_right * exp(log_1m_T));
    }
  }

  /* split Student-t log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: mean parameter
   *   sigma: positive scale parameter
   *   alpha: skewness parameter
   *   nu: degrees of freedom parameter (> 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real split_student_lccdf(real y, real mu, real sigma, real alpha, real nu) {
    real k = exp(alpha);
    real omega = sigma / split_s_k(alpha);
    real m = split_mode(mu, sigma, alpha, nu);
    real z = (y - m) / omega;
    real w_left  = 2.0 / (1.0 + square(k));
    real w_right = 2.0 * square(k) / (1.0 + square(k));
    
    if (z >= 0) {
      real u = z / k;
      // Guard at u = 0 ensures finite gradient for Stan's autodiff at the mode
      real log_1m_T = (u == 0) ? log(0.5) - 2 * u * exp(student_t_lpdf(0 | nu, 0, 1))
                                : student_t_lccdf(u | nu, 0, 1);
      return log(w_right) + log_1m_T;
    } else {
      real u = k * z;
      real log_T = (u == 0) ? log(0.5) + 2 * u * exp(student_t_lpdf(0 | nu, 0, 1))
                            : student_t_lcdf(u | nu, 0, 1);
      return log1m(w_left * exp(log_T));
    }
  }
