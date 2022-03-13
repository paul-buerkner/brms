  /* multi-logit transform
   * Args:
   *   y: simplex vector of length D
   *   ref: a single integer in 1:D indicating the reference category
   * Returns:
   *   an unbounded real vector of length D - 1
   */
   vector multi_logit(vector y, int ref) {
     vector[rows(y) - 1] x;
     for (i in 1:(ref - 1)) {
       x[i] = log(y[i]) - log(y[ref]);
     }
      for (i in (ref+1):rows(y)) {
       x[i - 1] = log(y[i]) - log(y[ref]);
     }
     return(x);
   }
  /* logistic-normal log-PDF
   * Args:
   *   y: simplex vector of response values (length D)
   *   mu: vector of means on the logit scale (length D-1)
   *   sigma: vector for standard deviations on the logit scale (length D-1)
   *   Lcor: Cholesky correlation matrix on the logit scale (dim D-1)
   *   ref: a single integer in 1:D indicating the reference category
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real logistic_normal_cholesky_cor_lpdf(vector y, vector mu, vector sigma,
                                          matrix Lcor, int ref) {
     int D = rows(y);
     vector[D - 1] x = multi_logit(y, ref);
     matrix[D - 1, D - 1] Lcov = diag_pre_multiply(sigma, Lcor);
     // multi-normal plus Jacobian adjustment of multivariate logit transform
     return multi_normal_cholesky_lpdf(x | mu, Lcov) - sum(log(y));
   }
