  /* compute the inverse of the cauchit link 
   * Args: 
   *   y: the real value to be transformed 
   * Returns: 
   *   a scalar in (0,1) 
   */ 
  real inv_cauchit(real y) { 
    real p; 
    p <- cauchy_cdf(y, 0, 1); 
    return p; 
  }