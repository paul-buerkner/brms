  /* Extended support beta log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   kappa: the exceedance parameter (u in Kosmidis & Zeileis, 2024, https://arxiv.org/abs/2409.07233)
   * Returns:
   *   a scalar to be added to the log posterior
   */

   real xbetax_lpdf(vector y, vector mu, vector phi, vector kappa) {
     int N = size(y);
     int N_zer = 0;
     int N_one = 0;
     int N_oth = 0;
     int i_zer = 0;
     int i_one = 0;
     int i_oth = 0;    
     for (i in 1:N) {
     	 N_zer += y[i] <= 0;
	 N_one += y[i] >= 1;
	 N_oth += (y[i] > 0) && (y[i] < 1);
     }     
     array[N_zer] int zer;
     array[N_one] int one;
     array[N_oth] int oth;
     for (i in 1:N) {
       if (y[i] <= 0) {
         i_zer += 1;
         zer[i_zer] = i;
       }
       else if (y[i] >= 1) {
         i_one += 1;
	 one[i_one] = i;
       }
       else {
	 i_oth += 1;
	 oth[i_oth] = i;
       }
     }
     real ll_zer = beta_proportion_lcdf(kappa[zer] ./ (1 + 2 * kappa[zer]) | mu[zer], phi[zer]);
     real ll_one = beta_proportion_lccdf((1 + kappa[one]) ./ (1 + 2 * kappa[one]) | mu[one], phi[one]);
     real ll_oth = beta_proportion_lpdf((y[oth] + kappa[oth]) ./ (1 + 2 * kappa[oth]) | mu[oth], phi[oth]) - sum(log(1 + 2 * kappa[oth]));
     real ll = ll_zer + ll_one + ll_oth;
     return ll;	
   }

   real xbetax_lpdf(vector y, vector mu, vector phi, real kappa) {
        int N = size(y);
     int N_zer = 0;
     int N_one = 0;
     int N_oth = 0;
     int i_zer = 0;
     int i_one = 0;
     int i_oth = 0;    
     for (i in 1:N) {
     	 N_zer += y[i] <= 0;
	 N_one += y[i] >= 1;
	 N_oth += (y[i] > 0) && (y[i] < 1);
     }     
     array[N_zer] int zer;
     array[N_one] int one;
     array[N_oth] int oth;
     for (i in 1:N) {
       if (y[i] <= 0) {
         i_zer += 1;
         zer[i_zer] = i;
       }
       else if (y[i] >= 1) {
         i_one += 1;
	 one[i_one] = i;
       }
       else {
	 i_oth += 1;
	 oth[i_oth] = i;
       }
     }
     real ll_zer = beta_proportion_lcdf(kappa / (1 + 2 * kappa) | mu[zer], phi[zer]);
     real ll_one = beta_proportion_lccdf((1 + kappa) / (1 + 2 * kappa) | mu[one], phi[one]);
     real ll_oth = beta_proportion_lpdf((y[oth] + kappa) / (1 + 2 * kappa) | mu[oth], phi[oth]) - N_oth * log(1 + 2 * kappa);
     real ll = ll_zer + ll_one + ll_oth;
     return ll;	
   }

   real xbetax_lpdf(vector y, vector mu, real phi, vector kappa) {
     int N = size(y);
     int N_zer = 0;
     int N_one = 0;
     int N_oth = 0;
     int i_zer = 0;
     int i_one = 0;
     int i_oth = 0;    
     for (i in 1:N) {
     	 N_zer += y[i] <= 0;
	 N_one += y[i] >= 1;
	 N_oth += (y[i] > 0) && (y[i] < 1);
     }     
     array[N_zer] int zer;
     array[N_one] int one;
     array[N_oth] int oth;
     for (i in 1:N) {
       if (y[i] <= 0) {
         i_zer += 1;
         zer[i_zer] = i;
       }
       else if (y[i] >= 1) {
         i_one += 1;
	 one[i_one] = i;
       }
       else {
	 i_oth += 1;
	 oth[i_oth] = i;
       }
     }
     real ll_zer = beta_proportion_lcdf(kappa[zer] ./ (1 + 2 * kappa[zer]) | mu[zer], phi);
     real ll_one = beta_proportion_lccdf((1 + kappa[one]) ./ (1 + 2 * kappa[one]) | mu[one], phi);
     real ll_oth = beta_proportion_lpdf((y[oth] + kappa[oth]) ./ (1 + 2 * kappa[oth]) | mu[oth], phi) - sum(log(1 + 2 * kappa[oth]));
     real ll = ll_zer + ll_one + ll_oth;
     return ll;	   
   }

  real xbetax_lpdf(vector y, vector mu, real phi, real kappa) {
     int N = size(y);
     int N_zer = 0;
     int N_one = 0;
     int N_oth = 0;
     int i_zer = 0;
     int i_one = 0;
     int i_oth = 0;    
     for (i in 1:N) {
     	 N_zer += y[i] <= 0;
	 N_one += y[i] >= 1;
	 N_oth += (y[i] > 0) && (y[i] < 1);
     }     
     array[N_zer] int zer;
     array[N_one] int one;
     array[N_oth] int oth;
     for (i in 1:N) {
       if (y[i] <= 0) {
         i_zer += 1;
         zer[i_zer] = i;
       }
       else if (y[i] >= 1) {
         i_one += 1;
	 one[i_one] = i;
       }
       else {
	 i_oth += 1;
	 oth[i_oth] = i;
       }
     }
     real ll_zer = beta_proportion_lcdf(kappa / (1 + 2 * kappa) | mu[zer], phi);
     real ll_one = beta_proportion_lccdf((1 + kappa) / (1 + 2 * kappa) | mu[one], phi);
     real ll_oth = beta_proportion_lpdf((y[oth] + kappa) / (1 + 2 * kappa) | mu[oth], phi) - N_oth * log(1 + 2 * kappa);
     real ll = ll_zer + ll_one + ll_oth;
     return ll;	
   }

   real xbetax_lpdf(real y, real mu, real phi, real kappa) {
     real denom = (1 + 2 * kappa);
     if (y <= 0)
       return beta_proportion_lcdf(kappa / denom | mu, phi);
     else if (y >= 1)
       return beta_proportion_lccdf((1 + kappa) / denom | mu, phi);
     else
       return beta_proportion_lpdf((y + kappa) / denom | mu, phi) - log(denom);
   }
