  /* Ordered beta regression log-PDF for a single response
   * Based on Kubinec (2023): https://doi.org/10.1017/pan.2022.20
   * Args:
   *   y: the response value in [0, 1]
   *   mu: latent mean parameter
   *   phi: precision parameter of the beta distribution
   *   thres: ordered thresholds (2 elements)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real ordbeta_logit_lpdf(real y, real mu, real phi, vector thres) {
    if (y == 0) {
      return log_inv_logit(thres[1] - mu);
    } else if (y == 1) {
      return log1m_inv_logit(thres[2] - mu);
    } else {
      real p_cont = inv_logit(thres[2] - mu) - inv_logit(thres[1] - mu);
      real beta_mu = inv_logit(mu);
      return log(p_cont) + beta_lpdf(y | beta_mu * phi, (1 - beta_mu) * phi);
    }
  }

  /* Ordered beta regression log-PDF for a single response (probit link)
   */
  real ordbeta_probit_lpdf(real y, real mu, real phi, vector thres) {
    if (y == 0) {
      return std_normal_lcdf(thres[1] - mu);
    } else if (y == 1) {
      return std_normal_lccdf(thres[2] - mu);
    } else {
      real p_cont = Phi(thres[2] - mu) - Phi(thres[1] - mu);
      real beta_mu = Phi(mu);
      return log(p_cont) + beta_lpdf(y | beta_mu * phi, (1 - beta_mu) * phi);
    }
  }

  /* Ordered beta regression log-PDF for a single response (cloglog link)
   */
  real ordbeta_cloglog_lpdf(real y, real mu, real phi, vector thres) {
    if (y == 0) {
      return log1m_exp(-exp(thres[1] - mu));
    } else if (y == 1) {
      return -exp(thres[2] - mu);
    } else {
      real p1 = inv_cloglog(thres[1] - mu);
      real p2 = inv_cloglog(thres[2] - mu);
      real p_cont = p2 - p1;
      real beta_mu = inv_cloglog(mu);
      return log(p_cont) + beta_lpdf(y | beta_mu * phi, (1 - beta_mu) * phi);
    }
  }

  /* Ordered beta regression log-PDF for a single response (cauchit link)
   */
  real ordbeta_cauchit_lpdf(real y, real mu, real phi, vector thres) {
    real p1 = cauchy_cdf(thres[1] - mu | 0, 1);
    real p2 = cauchy_cdf(thres[2] - mu | 0, 1);
    if (y == 0) {
      return log(p1);
    } else if (y == 1) {
      return log1m(p2);
    } else {
      real p_cont = p2 - p1;
      real beta_mu = cauchy_cdf(mu | 0, 1);
      return log(p_cont) + beta_lpdf(y | beta_mu * phi, (1 - beta_mu) * phi);
    }
  }

  // Vectorized version for logit link: y is vector, mu is vector, phi is scalar, thres is vector
  real ordbeta_logit_lpdf(vector y, vector mu, real phi, vector thres) {
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
    real ll_zer = sum(log_inv_logit(thres[1] - mu[zer]));
    real ll_one = sum(log1m_inv_logit(thres[2] - mu[one]));
    vector[N_oth] p_cont = inv_logit(thres[2] - mu[oth]) - inv_logit(thres[1] - mu[oth]);
    vector[N_oth] beta_mu = inv_logit(mu[oth]);
    real ll_oth = sum(log(p_cont)) + beta_lpdf(y[oth] | beta_mu * phi, (1 - beta_mu) * phi);
    return ll_zer + ll_one + ll_oth;
  }

  // Vectorized version for logit link: y is vector, mu is vector, phi is vector, thres is vector
  real ordbeta_logit_lpdf(vector y, vector mu, vector phi, vector thres) {
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
    real ll_zer = sum(log_inv_logit(thres[1] - mu[zer]));
    real ll_one = sum(log1m_inv_logit(thres[2] - mu[one]));
    vector[N_oth] p_cont = inv_logit(thres[2] - mu[oth]) - inv_logit(thres[1] - mu[oth]);
    vector[N_oth] beta_mu = inv_logit(mu[oth]);
    real ll_oth = sum(log(p_cont)) + beta_lpdf(y[oth] | beta_mu .* phi[oth], (1 - beta_mu) .* phi[oth]);
    return ll_zer + ll_one + ll_oth;
  }

  // Vectorized version for probit link
  real ordbeta_probit_lpdf(vector y, vector mu, real phi, vector thres) {
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
    real ll_zer = std_normal_lcdf(thres[1] - mu[zer] | );
    real ll_one = std_normal_lccdf(thres[2] - mu[one] | );
    vector[N_oth] p_cont = Phi(thres[2] - mu[oth]) - Phi(thres[1] - mu[oth]);
    vector[N_oth] beta_mu = Phi(mu[oth]);
    real ll_oth = sum(log(p_cont)) + beta_lpdf(y[oth] | beta_mu * phi, (1 - beta_mu) * phi);
    return ll_zer + ll_one + ll_oth;
  }

  // Vectorized version for probit link with vector phi
  real ordbeta_probit_lpdf(vector y, vector mu, vector phi, vector thres) {
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
    real ll_zer = std_normal_lcdf(thres[1] - mu[zer] | );
    real ll_one = std_normal_lccdf(thres[2] - mu[one] | );
    vector[N_oth] p_cont = Phi(thres[2] - mu[oth]) - Phi(thres[1] - mu[oth]);
    vector[N_oth] beta_mu = Phi(mu[oth]);
    real ll_oth = sum(log(p_cont)) + beta_lpdf(y[oth] | beta_mu .* phi[oth], (1 - beta_mu) .* phi[oth]);
    return ll_zer + ll_one + ll_oth;
  }
