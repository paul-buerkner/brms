  /* hurdle gamma log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   shape: shape parameter of gamma distribution 
   *   eta: linear predictor for gamma part 
   *   eta_hu: linear predictor for hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_gamma_lpdf(real y, real shape, real eta, real eta_hu) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | eta_hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | eta_hu) +  
              gamma_lpdf(y | shape, shape / exp(eta)); 
     } 
   }
