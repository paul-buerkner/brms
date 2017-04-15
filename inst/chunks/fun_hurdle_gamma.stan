  /* hurdle gamma log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   alpha: shape parameter of the gamma distribution 
   *   beta: rate parameter of the gamma distribution 
   *   hu: hurdle probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_gamma_lpdf(real y, real alpha, real beta, real hu) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              gamma_lpdf(y | alpha, beta); 
     } 
   }
  /* hurdle gamma log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   alpha: shape parameter of the gamma distribution 
   *   beta: rate parameter of the gamma distribution 
   *   hu: linear predictor for the hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_gamma_logit_lpdf(real y, real alpha, real beta, real hu) {
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | hu) +  
              gamma_lpdf(y | alpha, beta); 
     } 
   }
