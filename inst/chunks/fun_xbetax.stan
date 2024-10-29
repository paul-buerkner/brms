  /* Extended support beta log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   kappa: the exceedance parameter (u in Kosmidis & Zeileis, 2024, https://arxiv.org/abs/2409.07233)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   // real xbetax_lpdf(real y, real mu, real phi, real kappa) {
   //   real denom = (1 + 2 * kappa);
   //   if (y <= 0)
   //     return beta_proportion_lcdf(kappa / denom | mu, phi);
   //   else if (y >= 1)
   //     return beta_proportion_lccdf((1 + kappa) / denom | mu, phi);
   //   else
   //     return beta_proportion_lpdf((y + kappa) / denom | mu, phi) - log(denom);
   // }

   real xbetax_lpdf(vector y, vector mu, vector phi, vector kappa) {
     int N;
     N = size(y);
     vector[N] out;
     vector[N] denom = (1 + 2 * kappa);
     real ll = 0;
     for (i in 1:N) {
       if (y[i] <= 0) {
         ll = ll + beta_proportion_lcdf(kappa[i] / denom[i] | mu[i], phi[i]);
       }
       else {
         if (y[i] >= 1) {
           ll = ll + beta_proportion_lccdf((1 + kappa[i]) / denom[i] | mu[i], phi[i]);
	 }
         else {
	   ll = ll + beta_proportion_lpdf((y[i] + kappa[i]) / denom[i] | mu[i], phi[i]) - log(denom[i]);
	 }
       } 
     }
     return ll;	
   }

   real xbetax_rng(real mu, real phi, real kappa) {
     real z;
     z = (1 + 2 * kappa) * beta_proportion_rng(mu, phi) - kappa;
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
