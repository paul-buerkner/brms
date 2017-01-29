  /* compute the logit_m1_to_half link 
   * Args: 
   *   p: a scalar in (-1, 0.5)
   * Returns: 
   *   a scalar in (-Inf, Inf)
   */ 
   real logit_m1_to_half(real p) { 
     return logit((p + 1) / 1.5); 
   }
  /* compute the inverse of the logit_m1_to_half link 
   * Args: 
   *   y: a scalar in (-Inf, Inf)
   * Returns: 
   *   a scalar in (-1, 0.5)
   */ 
   real inv_logit_m1_to_half(real y) { 
     return inv_logit(y) * 1.5 - 1;
   }
