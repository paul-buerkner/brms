/* shifted log-Student-t log-PDF for a single response
 * Args:
 *   y: response value (must be > ndt)
 *   mu: location parameter on the log scale
 *   sigma: positive scale parameter on the log scale
 *   nu: degrees of freedom
 *   ndt: non-decision time (shift)
 * Returns:
 *   a scalar to be added to the log posterior
 */
real shifted_logstudent_lpdf(real y, real mu, real sigma, real nu, real ndt) {
  if (y <= ndt) {
    return negative_infinity();
  }
  real log_diff = log(y - ndt);
  return student_t_lpdf(log_diff | nu, mu, sigma) - log_diff;
}

// shifted log-Student-t log-CDF
real shifted_logstudent_lcdf(real y, real mu, real sigma, real nu, real ndt) {
  if (y <= ndt) {
    return negative_infinity();
  }
  return student_t_lcdf(log(y - ndt) | nu, mu, sigma);
}

// shifted log-Student-t log-CCDF
real shifted_logstudent_lccdf(real y, real mu, real sigma, real nu, real ndt) {
  if (y <= ndt) {
    return 0.0;
  }
  return student_t_lccdf(log(y - ndt) | nu, mu, sigma);
}
