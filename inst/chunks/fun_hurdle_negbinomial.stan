  /* hurdle negative binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   eta_hu: linear predictor for hurdle part 
   *   shape: shape parameter of negative binomial distribution 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_neg_binomial_2_lpmf(int y, real eta, real eta_hu,  
                                   real shape) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | eta_hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | eta_hu) +  
              neg_binomial_2_log_lpmf(y | eta, shape) - 
              log(1 - (shape / (exp(eta) + shape))^shape); 
     } 
   }
