  /* Extended support beta log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   u: the exceedance parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real xbetax_lpdf(real y, real mu, real phi, real u) {
     real denom = (1 + 2 * u);
     if (y <= 0)
       return beta_proportion_lcdf(u / denom | mu, phi);
     else if (y >= 1)
       return beta_proportion_lccdf((1 + u) / denom | mu, phi);
     else
       return beta_proportion_lpdf((y + u) / denom | mu, phi) - log(denom);
   }

   real xbetax_rng(real mu, real phi, real u) {
     real z;
     z = (1 + 2 * u) * beta_proportion_rng(mu, phi) - u;
     if (z < 0) {
       return 0;
     } else {
       if (z > 1) {
         return 1;
       } else {
         return z;
       }
     }
   }
