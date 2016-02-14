  /* zero-inflated negative binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   eta_zi: linear predictor for zero-inflation part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_neg_binomial_2_log(int y, real eta, real eta_zi, 
                                         real shape) { 
     if (y == 0) { 
       return log_sum_exp(bernoulli_logit_log(1, eta_zi), 
                          bernoulli_logit_log(0, eta_zi) + 
                          neg_binomial_2_log_log(0, eta, shape)); 
     } else { 
       return bernoulli_logit_log(0, eta_zi) +  
              neg_binomial_2_log_log(y, eta, shape); 
     } 
   } 
