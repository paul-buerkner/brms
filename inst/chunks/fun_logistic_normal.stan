  /* multi-logit transform
   * Args: 
   *   y: a simplex vector
   * Returns:  
   *   an unbounded real vector of length size(y) - 1
   */ 
   vector multi_logit(vector y, int ref) {
     vector[size(y) - 1] x;
     for (i in 1:(ref - 1)) {
       x[i] = log(y[i] / y[ref]);
     }
      for (i in (ref+1):size(y)) {
       x[i - 1] = log(y[i] / y[ref]); 
     }
     return(x);
   } 
  /* logistic-normal log-PDF
   * Args: 
   *   y: vector of real response values
   *   mu: vector of means of the logit scale
   *   sigma: vector for standard deviations on the logit scale
   *   Lcor: correlation matrix on the logit scale
   *   ref: index of the mandatory reference category
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real logistic_normal_cholesky_cor_lpdf(vector y, vector mu, vector sigma, 
                                          matrix Lcor, int ref) {
     int D = size(y);
     vector[D - 1] x = multi_logit(y, ref);
     matrix[D - 1, D - 1] Lcov = diag_pre_multiply(sigma, Lcor);
     // multi-normal plus Jacobian adjustment of multivariate logit transform
     return multi_normal_cholesky_lpdf(x | mu, Lcov) - sum(log(y));
   }
