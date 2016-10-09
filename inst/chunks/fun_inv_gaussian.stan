  /* inverse Gaussian log-PDF for a single response (for data only) 
   * Copyright Stan Development Team 2015 
   * Args: 
   *   y: the response value 
   *   mu: positive mean parameter 
   *   shape: positive shape parameter 
   *   log_y: precomputed log(y) 
   *   sqrt_y: precomputed sqrt(y) 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real inv_gaussian_lpdf(real y, real mu, real shape,  
                          real log_y, real sqrt_y) { 
     return 0.5 * log(shape / (2 * pi())) -  
            1.5 * log_y - 
            0.5 * shape * square((y - mu) / (mu * sqrt_y)); 
   }
  /* vectorized inverse Gaussian log-PDF (for data only) 
   * Copyright Stan Development Team 2015 
   * Args: 
   *   y: response vector 
   *   mu: positive mean parameter vector 
   *   shape: positive shape parameter 
   *   sum_log_y: precomputed sum of log(y) 
   *   sqrt_y: precomputed sqrt(y) 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real inv_gaussian_vector_lpdf(vector y, vector mu, real shape,  
                                 real sum_log_y, vector sqrt_y) { 
     return 0.5 * rows(y) * log(shape / (2 * pi())) -  
            1.5 * sum_log_y - 
            0.5 * shape * dot_self((y - mu) ./ (mu .* sqrt_y)); 
   }
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
   real inv_gaussian_lcdf(real y, real mu, real shape,  
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
   real inv_gaussian_lccdf(real y, real mu, real shape, 
                           real log_y, real sqrt_y) { 
     return log(1 - Phi(sqrt(shape) / sqrt_y * (y / mu - 1)) - 
                  exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt_y * (y / mu + 1)));
   }
