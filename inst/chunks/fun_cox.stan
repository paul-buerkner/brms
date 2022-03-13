  /* distribution functions of the Cox proportional hazards model
   * parameterize hazard(t) = baseline(t) * mu
   * so that higher values of 'mu' imply lower survival times
   * Args:
   *   y: the response value; currently ignored as the relevant
   *     information is passed via 'bhaz' and 'cbhaz'
   *   mu: positive location parameter
   *   bhaz: baseline hazard
   *   cbhaz: cumulative baseline hazard
   */
  real cox_lhaz(real y, real mu, real bhaz, real cbhaz) {
    return log(bhaz) + log(mu);
  }
  real cox_lccdf(real y, real mu, real bhaz, real cbhaz) {
    // equivalent to the log survival function
    return - cbhaz * mu;
  }
  real cox_lcdf(real y, real mu, real bhaz, real cbhaz) {
    return log1m_exp(cox_lccdf(y | mu, bhaz, cbhaz));
  }
  real cox_lpdf(real y, real mu, real bhaz, real cbhaz) {
    return cox_lhaz(y, mu, bhaz, cbhaz) + cox_lccdf(y | mu, bhaz, cbhaz);
  }
  // Distribution functions of the Cox model in log parameterization
  real cox_log_lhaz(real y, real log_mu, real bhaz, real cbhaz) {
    return log(bhaz) + log_mu;
  }
  real cox_log_lccdf(real y, real log_mu, real bhaz, real cbhaz) {
    return - cbhaz * exp(log_mu);
  }
  real cox_log_lcdf(real y, real log_mu, real bhaz, real cbhaz) {
    return log1m_exp(cox_log_lccdf(y | log_mu, bhaz, cbhaz));
  }
  real cox_log_lpdf(real y, real log_mu, real bhaz, real cbhaz) {
    return cox_log_lhaz(y, log_mu, bhaz, cbhaz) +
           cox_log_lccdf(y | log_mu, bhaz, cbhaz);
  }
