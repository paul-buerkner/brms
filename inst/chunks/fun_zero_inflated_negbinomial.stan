  /* zero-inflated negative binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   zi: zero-inflation probability
   *   shape: shape parameter of negative binomial distribution
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_neg_binomial_lpmf(int y, real eta, real zi, 
                                          real shape) { 
     if (y == 0) { 
       return log_sum_exp(bernoulli_lpmf(1 | zi), 
                          bernoulli_lpmf(0 | zi) + 
                          neg_binomial_2_log_lpmf(0 | eta, shape)); 
     } else { 
       return bernoulli_lpmf(0 | zi) +  
              neg_binomial_2_log_lpmf(y | eta, shape); 
     } 
   } 
  /* zero-inflated negative binomial log-PDF of a single response 
   * logit parameterization of the zero-inflation part
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   zi: linear predictor for zero-inflation part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_neg_binomial_logit_lpmf(int y, real eta, 
                                              real zi, real shape) { 
     if (y == 0) { 
       return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                          bernoulli_logit_lpmf(0 | zi) + 
                          neg_binomial_2_log_lpmf(0 | eta, shape)); 
     } else { 
       return bernoulli_logit_lpmf(0 | zi) +  
              neg_binomial_2_log_lpmf(y | eta, shape); 
     } 
   }
