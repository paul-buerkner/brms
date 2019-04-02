 /* discrete Weibull log-PMF for a single response
  * Args: 
  *   y: the response value 
  *   mu: location parameter on the unit interval
  *   shape: positive shape parameter
  * Returns:  
  *   a scalar to be added to the log posterior 
  */ 
  real discrete_weibull_lpmf(int y, real mu, real shape) {
    return log(mu^y^shape - mu^(y+1)^shape);
  }
  // discrete Weibull log-CDF for a single response
  real discrete_weibull_lcdf(int y, real mu, real shape) {
    return log(1 - mu^(y + 1)^shape);
  }
  // discrete Weibull log-CCDF for a single response
  real discrete_weibull_lccdf(int y, real mu, real shape) {
    return log(mu) * (y + 1)^shape;
  }
