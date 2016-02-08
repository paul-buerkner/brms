  /* hurdle poisson log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for poisson part 
   *   eta_hu: linear predictor for hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_poisson_log(int y, real eta, real eta_hu) { 
     if (y == 0) { 
       return bernoulli_logit_log(1, eta_hu); 
     } else { 
       return bernoulli_logit_log(0, eta_hu) +  
              poisson_log_log(y, eta) - 
              log(1 - exp(-exp(eta))); 
     } 
   }
  