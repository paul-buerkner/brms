  /* exponentially modified Gaussian log-PDF for a single response
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the gaussian component
   *   sigma: SD parameter of the gaussian component
   *   beta: scale parameter of the exponential component
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real exgaussian_lpdf(real y, real mu, real sigma, real beta) { 
     real z;
     z = y - mu - sigma^2 / beta;
     if (beta > 0.05 * sigma) {
       return -log(beta) - (z + sigma^2 / (2 * beta)) / beta + log(Phi(z / sigma));
     } else {
       return normal_lpdf(y | mu, sigma);
     }
   }
  /* exponentially modified Gaussian log-CDF for a single quantile
   * Args: 
   *   y: a quantile
   *   mu: mean parameter of the gaussian component
   *   sigma: SD parameter of the gaussian component
   *   beta: scale parameter of the exponential component
   * Returns:  
   *   log(P(Y <= y))
   */ 
   real exgaussian_lcdf(real y, real mu, real sigma, real beta) { 
     real z;
     z = y - mu - sigma^2 / beta;
     if (beta > 0.05 * sigma) {
       return log(Phi((y - mu) / sigma) - Phi(z / sigma) * 
         exp(((mu + sigma^2 / beta)^2 - mu^2 - 2 * y * sigma^2 / beta) / 
             (2 * sigma^2)));
     } else {
       return normal_lcdf(y | mu, sigma);
     }
   }
  /* exponentially modified Gaussian log-CCDF for a single quantile
   * Args: 
   *   y: a quantile
   *   mu: mean parameter of the gaussian component
   *   sigma: SD parameter of the gaussian component
   *   beta: scale parameter of the exponential component
   * Returns:  
   *   log(P(Y > y))
   */ 
   real exgaussian_lccdf(real y, real mu, real sigma, real beta) { 
     real z;
     z = y - mu - sigma^2 / beta;
     if (beta > 0.05 * sigma) {
       return log(1 - (Phi((y - mu) / sigma) - Phi(z / sigma) * 
         exp(((mu + sigma^2 / beta)^2 - mu^2 - 2 * y * sigma^2 / beta) / 
             (2 * sigma^2))));
     } else {
       return normal_lccdf(y | mu, sigma);
     }
   }
