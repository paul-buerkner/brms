  /* Ordered beta regression log-PDF for a single response
   * Based on Kubinec (2023): https://doi.org/10.1017/pan.2022.20
   * Args:
   *   y: the response value in [0, 1]
   *   mu: mean parameter (on response scale, i.e., in (0, 1))
   *   phi: precision parameter of the beta distribution
   *   zoi: first cutpoint/threshold for boundary at 0 (latent scale)
   *   kappa: positive distance to second cutpoint, so coi = zoi + kappa
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real ordbeta_lpdf(real y, real mu, real phi, real zoi, real kappa) {
    // Compute coi (second cutpoint) from zoi and kappa
    real coi = zoi + kappa;
    // Transform mu to latent scale for threshold comparison
    real mu_latent = logit(mu);
    // Compute probabilities for each component
    real pr_zero = inv_logit(zoi - mu_latent);
    real pr_one = 1 - inv_logit(coi - mu_latent);
    real pr_cont = inv_logit(coi - mu_latent) - inv_logit(zoi - mu_latent);
    if (y == 0) {
      return log(pr_zero);
    } else if (y == 1) {
      return log(pr_one);
    } else {
      return log(pr_cont) + beta_lpdf(y | mu * phi, (1 - mu) * phi);
    }
  }

  // Vectorized version: y is vector, mu is vector, phi is scalar
  real ordbeta_lpdf(vector y, vector mu, real phi, real zoi, real kappa) {
    // Compute coi (second cutpoint) from zoi and kappa
    real coi = zoi + kappa;
    int N = size(y);
    int N_zer = 0;
    int N_one = 0;
    int N_oth = 0;
    int i_zer = 0;
    int i_one = 0;
    int i_oth = 0;

    for (i in 1:N) {
      N_zer += y[i] == 0;
      N_one += y[i] == 1;
      N_oth += (y[i] > 0) && (y[i] < 1);
    }
    array[N_zer] int zer;
    array[N_one] int one;
    array[N_oth] int oth;
    for (i in 1:N) {
      if (y[i] == 0) {
        i_zer += 1;
        zer[i_zer] = i;
      } else if (y[i] == 1) {
        i_one += 1;
        one[i_one] = i;
      } else {
        i_oth += 1;
        oth[i_oth] = i;
      }
    }
    // Transform mu to latent scale
    vector[N] mu_latent = logit(mu);
    real ll_zer = sum(log_inv_logit(zoi - mu_latent[zer]));
    real ll_one = sum(log1m_inv_logit(coi - mu_latent[one]));
    vector[N_oth] pr_cont = inv_logit(coi - mu_latent[oth]) - inv_logit(zoi - mu_latent[oth]);
    real ll_oth = sum(log(pr_cont)) + beta_lpdf(y[oth] | mu[oth] * phi, (1 - mu[oth]) * phi);
    return ll_zer + ll_one + ll_oth;
  }

  // Vectorized version: y is vector, mu is vector, phi is vector
  real ordbeta_lpdf(vector y, vector mu, vector phi, real zoi, real kappa) {
    // Compute coi (second cutpoint) from zoi and kappa
    real coi = zoi + kappa;
    int N = size(y);
    int N_zer = 0;
    int N_one = 0;
    int N_oth = 0;
    int i_zer = 0;
    int i_one = 0;
    int i_oth = 0;

    for (i in 1:N) {
      N_zer += y[i] == 0;
      N_one += y[i] == 1;
      N_oth += (y[i] > 0) && (y[i] < 1);
    }
    array[N_zer] int zer;
    array[N_one] int one;
    array[N_oth] int oth;
    for (i in 1:N) {
      if (y[i] == 0) {
        i_zer += 1;
        zer[i_zer] = i;
      } else if (y[i] == 1) {
        i_one += 1;
        one[i_one] = i;
      } else {
        i_oth += 1;
        oth[i_oth] = i;
      }
    }
    // Transform mu to latent scale
    vector[N] mu_latent = logit(mu);
    real ll_zer = sum(log_inv_logit(zoi - mu_latent[zer]));
    real ll_one = sum(log1m_inv_logit(coi - mu_latent[one]));
    vector[N_oth] pr_cont = inv_logit(coi - mu_latent[oth]) - inv_logit(zoi - mu_latent[oth]);
    real ll_oth = sum(log(pr_cont)) + beta_lpdf(y[oth] | mu[oth] .* phi[oth], (1 - mu[oth]) .* phi[oth]);
    return ll_zer + ll_one + ll_oth;
  }
