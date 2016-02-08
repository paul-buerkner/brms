  /* inverse Gaussian log-CDF for a single quantile 
   * Args: 
   *   y: a quantile 
   *   mu: positive mean parameter 
   *   shape: positive shape parameter 
   *   log_y: ignored (cdf and pdf should have the same args) 
   *   sqrt_y: precomputed sqrt(y) 
   * Returns: 
   *   log(P(Y <= y)) 
   */ 
  real inv_gaussian_cdf_log(real y, real mu, real shape,  
                            real log_y, real sqrt_y) { 
    return log(Phi(sqrt(shape) / sqrt_y * (y / mu - 1)) + 
               exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt_y * (y / mu + 1))); 
  } 
  /* inverse Gaussian log-CCDF for a single quantile 
   * Args: 
   *   y: a quantile 
   *   mu: positive mean parameter 
   *   shape: positive shape parameter 
   *   log_y: ignored (ccdf and pdf should have the same args) 
   *   sqrt_y: precomputed sqrt(y) 
   * Returns: 
   *   log(P(Y > y)) 
   */ 
  real inv_gaussian_ccdf_log(real y, real mu, real shape, 
                             real log_y, real sqrt_y) { 
    return log(1 - Phi(sqrt(shape) / sqrt_y * (y / mu - 1)) - 
               exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt_y * (y / mu + 1))); 
  }