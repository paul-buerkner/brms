  /* First absolute moment factor b_nu for Azzalini skew-t
   * Args:
   *   nu: degrees of freedom parameter (> 1)
   */
  real st_b_nu(real nu) {
    return sqrt(nu) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2)) / sqrt(pi());
  }

  /* Azzalini-Capitanio skew Student-t log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter
   *   sigma: positive scale parameter (asymptotic SD as nu -> Inf)
   *   alpha: skewness parameter
   *   nu: degrees of freedom parameter (> 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real skew_student_lpdf(real y, real mu, real sigma, real alpha, real nu) {
    real delta = alpha / sqrt(1 + square(alpha));
    real omega = sigma / sqrt(1 - 2 * square(delta) / pi());
    real z = (y - mu) / omega + delta * st_b_nu(nu);
    real x = alpha * z * sqrt((nu + 1) / (nu + square(z)));
    // Stan's student_t_lcdf has a non-finite gradient at exactly x = 0 (e.g. alpha = 0 with init = 0).
    // There, use log(1/2) + 2 * x * exp(student_t_lpdf(0 | nu + 1, 0, 1)): same value, correct derivatives.
    real log_T = (x == 0) ? log(0.5) + 2 * x * exp(student_t_lpdf(0 | nu + 1, 0, 1))
                          : student_t_lcdf(x | nu + 1, 0, 1);
    return log(2.0) - log(omega) + student_t_lpdf(z | nu, 0, 1) + log_T;
  }
