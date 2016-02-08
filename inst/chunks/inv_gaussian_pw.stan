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
   real inv_gaussian_log(real y, real mu, real shape,  
                         real log_y, real sqrt_y) { 
     return 0.5 * log(shape / (2 * pi())) -  
            1.5 * log_y - 
            0.5 * shape * square((y - mu) / (mu * sqrt_y)); 
   }