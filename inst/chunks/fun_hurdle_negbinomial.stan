  /* hurdle negative binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   hu: hurdle probability
   *   shape: shape parameter of negative binomial distribution
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_neg_binomial_lpmf(int y, real eta, 
                                 real hu, real shape) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              neg_binomial_2_log_lpmf(y | eta, shape) - 
              log(1 - (shape / (exp(eta) + shape))^shape); 
     } 
   }
  /* hurdle negative binomial log-PDF of a single response 
   * logit parameterization for the hurdle part
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for negative binomial part 
   *   hu: linear predictor of hurdle part 
   *   shape: shape parameter of negative binomial distribution 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_neg_binomial_logit_lpmf(int y, real eta, 
                                       real hu, real shape) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | hu) +  
              neg_binomial_2_log_lpmf(y | eta, shape) - 
              log(1 - (shape / (exp(eta) + shape))^shape); 
     } 
   }
