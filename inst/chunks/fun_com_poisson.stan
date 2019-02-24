 /* log normalizing constant of the COM Poisson distribution
  * Args:
  *   log_mu: log location parameter
  *   shape: positive shape parameter
  */
  real log_Z_com_poisson(real log_mu, real nu) {
    real log_Z;
    real lfac;
    real term;
    real k;
    int M;
    real log_thres;
    if (nu == 1) {
      return exp(log_mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite")
    }
    M = 10000;
    log_thres = log(1e-16);
    // check if the Mth term of the series is small enough
    if (nu * (M * log_mu - lgamma(M + 1)) > log_thres) {
      reject("nu is too close to zero.")
    }
    log_Z = log1p_exp(nu * log_mu);  // first 2 terms of the series
    lfac = 0;
    term = 0;
    k = 2;
    while (term > log_thres) { 
      lfac += log(k);
      term = nu * (k * log_mu - lfac);
      log_Z = log_sum_exp(log_Z, term);
      k += 1;
    }
    return log_Z;
  }
 /* COM Poisson log-PMF for a single response (log parameterization)
  * Args: 
  *   y: the response value 
  *   log_mu: log location parameter
  *   shape: positive shape parameter
  */ 
  real com_poisson_log_lpmf(int y, real log_mu, real nu) {
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return nu * (y * log_mu - lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu);
  }
  // COM Poisson log-PMF for a single response
  real com_poisson_lpmf(int y, real mu, real nu) {
    if (nu == 1) return poisson_lpmf(y | mu);
    return com_poisson_log_lpmf(y | log(mu), nu);
  }
  // COM Poisson log-CDF for a single response
  real com_poisson_lcdf(int y, real mu, real nu) {
    int M;
    real log_thres;
    real log_mu;
    real log_num;  // log numerator
    real log_Z;  // log denominator
    if (nu == 1) {
      return poisson_lcdf(y | mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite")
    }
    M = 10000;
    if (y > M) {
      reject("cannot handle y > 10000")
    }
    log_thres = log(1e-16);
    log_mu = log(mu);
    if (nu * (y * log_mu - lgamma(y + 1)) <= log_thres) {
      // y is large enough for the CDF to be very close to 1;
      return 0;
    }
    if (y == 0) {
      return -log_Z_com_poisson(log_mu, nu); 
    } 
    if (y == 1) {
      return log1p_exp(nu * log_mu) - log_Z_com_poisson(log_mu, nu);
    }
    // check if the Mth term of the series is small enough
    if (nu * (M * log_mu - lgamma(M + 1)) > log_thres) {
      reject("nu is too close to zero.")
    }
    log_Z = log1p_exp(nu * log_mu);  // first 2 terms of the series
    {
      real lfac = 0;
      real term = 0;
      real k = 2;
      while (k <= y) {
        lfac += log(k);
        term = nu * (k * log_mu - lfac);
        log_Z = log_sum_exp(log_Z, term);
        k += 1;
      }
      // log_num consists of the first y + 1 terms of log_Z
      log_num = log_Z; 
      while (term > log_thres) { 
        lfac += log(k);
        term = nu * (k * log_mu - lfac);
        log_Z = log_sum_exp(log_Z, term);
        k += 1;
      }
    }
    return log_num - log_Z;
  }
  // COM Poisson log-CCDF for a single response
  real com_poisson_lccdf(int y, real mu, real nu) {
    return log1m_exp(com_poisson_lcdf(y | mu, nu));   
  }
