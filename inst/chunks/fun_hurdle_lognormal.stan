  /* hurdle lognormal log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for lognormal part 
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_lpdf(real y, real eta, real hu, real sigma) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              lognormal_lpdf(y | eta, sigma); 
     } 
   }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for lognormal part 
   *   eta_hu: linear predictor for hurdle part 
   *   sigma: sd parameter of the lognormal distribution
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_logit_lpdf(real y, real eta, real eta_hu, real sigma) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | eta_hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | eta_hu) +  
              lognormal_lpdf(y | eta, sigma); 
     } 
   }
