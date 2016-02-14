  /* zero-inflated poisson log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for poisson part 
   *   eta_zi: linear predictor for zero-inflation part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_poisson_log(int y, real eta, real eta_zi) { 
     if (y == 0) { 
       return log_sum_exp(bernoulli_logit_log(1, eta_zi), 
                          bernoulli_logit_log(0, eta_zi) + 
                          poisson_log_log(0, eta)); 
     } else { 
       return bernoulli_logit_log(0, eta_zi) +  
              poisson_log_log(y, eta); 
     } 
   }
