  /* zero-inflated beta log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_beta_lpdf(real y, real mu, real phi, real zi) {
     vector[2] shape; 
     shape[1] = mu * phi; 
     shape[2] = (1 - mu) * phi; 
     if (y == 0) { 
       return bernoulli_lpmf(1 | zi); 
     } else { 
       return bernoulli_lpmf(0 | zi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
  /* zero-inflated beta log-PDF of a single response 
   * logit parameterization of the zero-inflation part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_beta_logit_lpdf(real y, real mu, real phi, real zi) {
     vector[2] shape;
     shape[1] = mu * phi; 
     shape[2] = (1 - mu) * phi; 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | zi); 
     } else { 
       return bernoulli_logit_lpmf(0 | zi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
