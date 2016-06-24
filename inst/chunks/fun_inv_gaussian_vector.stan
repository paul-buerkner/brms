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
