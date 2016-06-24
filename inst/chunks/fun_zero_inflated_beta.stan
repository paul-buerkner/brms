  /* zero-inflated beta log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for beta part 
   *   eta_zi: linear predictor for zero-inflated part 
   *   phi: precision parameter 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_beta_lpdf(real y, real eta, real eta_zi, 
                                real phi) { 
     real inv_logit_eta; 
     vector[2] shape; 
     inv_logit_eta = inv_logit(eta); 
     shape[1] = inv_logit_eta * phi; 
     shape[2] = (1 - inv_logit_eta) * phi; 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | eta_zi); 
     } else { 
       return bernoulli_logit_lpmf(0 | eta_zi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
