  /* hurdle gamma log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   shape: shape parameter of gamma distribution 
   *   eta: linear predictor for gamma part 
   *   hu: hurdle probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_gamma_lpdf(real y, real shape, real eta, real hu) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              gamma_lpdf(y | shape, shape / exp(eta)); 
     } 
   }
  /* hurdle gamma log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   shape: shape parameter of gamma distribution 
   *   eta: linear predictor for gamma part 
   *   hu: linear predictor for hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_gamma_logit_lpdf(real y, real shape, 
                                real eta, real hu) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | hu) +  
              gamma_lpdf(y | shape, shape / exp(eta)); 
     } 
   }
