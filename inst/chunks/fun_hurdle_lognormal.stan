  /* hurdle lognormal log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the lognormal distribution 
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_lpdf(real y, real mu, real sigma, real hu) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              lognormal_lpdf(y | mu, sigma); 
     } 
   }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the lognormal distribution 
   *   sigma: sd parameter of the lognormal distribution
   *   hu: linear predictor for the hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_logit_lpdf(real y, real mu, real sigma, real hu) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | hu) +  
              lognormal_lpdf(y | mu, sigma); 
     } 
   }
