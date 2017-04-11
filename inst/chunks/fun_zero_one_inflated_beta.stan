  /* zero-one-inflated beta log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   eta: linear predictor for the beta part 
   *   zoi: zero-one-inflation probability
   *   coi: conditional one-inflation probability
   *   phi: precision parameter 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_one_inflated_beta_lpdf(real y, real eta, real zoi,
                                    real coi, real phi) {
     real inv_logit_eta;
     vector[2] shape;
     inv_logit_eta = inv_logit(eta); 
     shape[1] = inv_logit_eta * phi; 
     shape[2] = (1 - inv_logit_eta) * phi; 
     if (y == 0) { 
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(0 | coi); 
     } else if (y == 1) {
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(1 | coi);
     } else { 
       return bernoulli_lpmf(0 | zoi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
