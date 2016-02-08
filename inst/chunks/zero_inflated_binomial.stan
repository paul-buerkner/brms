  /* zero-inflated binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for binomial part 
   *   eta_zi: linear predictor for zero-inflation part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_binomial_log(int y, int trials, real eta, 
                                   real eta_zi) { 
     if (y == 0) { 
       return log_sum_exp(bernoulli_logit_log(1, eta_zi), 
                          bernoulli_logit_log(0, eta_zi) + 
                          binomial_logit_log(0, trials, eta)); 
     } else { 
       return bernoulli_logit_log(0, eta_zi) +  
              binomial_logit_log(y, trials, eta); 
     } 
   }
  