  /* Ordered beta regression log-PDF of a single response
   * Based on Kubinec (2023): https://doi.org/10.1017/pan.2022.20
   * Args:
   *   y: the response value in [0, 1]
   *   mu: linear predictor (on identity scale)
   *   phi: precision parameter of the beta distribution
   *   cutzero: first cutpoint parameter
   *   cutone: second cutpoint parameter (parameterized as log-offset from cutzero)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real ordbeta_lpdf(real y, real mu, real phi, real cutzero, real cutone) {
    real thresh1 = cutzero;
    real thresh2 = cutzero + exp(cutone);
    if (y == 0) {
      return log1m_inv_logit(mu - thresh1);
    } else if (y == 1) {
      return log_inv_logit(mu - thresh2);
    } else {
      return log_diff_exp(log_inv_logit(mu - thresh1), log_inv_logit(mu - thresh2)) +
             beta_lpdf(y | inv_logit(mu) * phi, (1 - inv_logit(mu)) * phi);
    }
  }

  // vectorized version: all parameters are vectors
  real ordbeta_lpdf(vector y, vector mu, vector phi, vector cutzero, vector cutone) {
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
    vector[N_zer] thresh1_zer = cutzero[zer];
    vector[N_one] thresh1_one = cutzero[one];
    vector[N_one] thresh2_one = cutzero[one] + exp(cutone[one]);
    vector[N_oth] thresh1_oth = cutzero[oth];
    vector[N_oth] thresh2_oth = cutzero[oth] + exp(cutone[oth]);
    real ll_zer = sum(log1m_inv_logit(mu[zer] - thresh1_zer));
    real ll_one = sum(log_inv_logit(mu[one] - thresh2_one));
    real ll_oth = sum(log_diff_exp(log_inv_logit(mu[oth] - thresh1_oth),
                                   log_inv_logit(mu[oth] - thresh2_oth))) +
                  beta_lpdf(y[oth] | inv_logit(mu[oth]) .* phi[oth],
                            (1 - inv_logit(mu[oth])) .* phi[oth]);
    return ll_zer + ll_one + ll_oth;
  }

  // vectorized version: mu is vector, phi/cutzero/cutone are scalars
  real ordbeta_lpdf(vector y, vector mu, real phi, real cutzero, real cutone) {
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
    real thresh1 = cutzero;
    real thresh2 = cutzero + exp(cutone);
    real ll_zer = sum(log1m_inv_logit(mu[zer] - thresh1));
    real ll_one = sum(log_inv_logit(mu[one] - thresh2));
    real ll_oth = sum(log_diff_exp(log_inv_logit(mu[oth] - thresh1),
                                   log_inv_logit(mu[oth] - thresh2))) +
                  beta_lpdf(y[oth] | inv_logit(mu[oth]) * phi,
                            (1 - inv_logit(mu[oth])) * phi);
    return ll_zer + ll_one + ll_oth;
  }

  // vectorized version: mu/phi are vectors, cutzero/cutone are scalars
  real ordbeta_lpdf(vector y, vector mu, vector phi, real cutzero, real cutone) {
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
    real thresh1 = cutzero;
    real thresh2 = cutzero + exp(cutone);
    real ll_zer = sum(log1m_inv_logit(mu[zer] - thresh1));
    real ll_one = sum(log_inv_logit(mu[one] - thresh2));
    real ll_oth = sum(log_diff_exp(log_inv_logit(mu[oth] - thresh1),
                                   log_inv_logit(mu[oth] - thresh2))) +
                  beta_lpdf(y[oth] | inv_logit(mu[oth]) .* phi[oth],
                            (1 - inv_logit(mu[oth])) .* phi[oth]);
    return ll_zer + ll_one + ll_oth;
  }

  // vectorized version: mu/cutzero/cutone are vectors, phi is scalar
  real ordbeta_lpdf(vector y, vector mu, real phi, vector cutzero, vector cutone) {
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
    vector[N_zer] thresh1_zer = cutzero[zer];
    vector[N_one] thresh1_one = cutzero[one];
    vector[N_one] thresh2_one = cutzero[one] + exp(cutone[one]);
    vector[N_oth] thresh1_oth = cutzero[oth];
    vector[N_oth] thresh2_oth = cutzero[oth] + exp(cutone[oth]);
    real ll_zer = sum(log1m_inv_logit(mu[zer] - thresh1_zer));
    real ll_one = sum(log_inv_logit(mu[one] - thresh2_one));
    real ll_oth = sum(log_diff_exp(log_inv_logit(mu[oth] - thresh1_oth),
                                   log_inv_logit(mu[oth] - thresh2_oth))) +
                  beta_lpdf(y[oth] | inv_logit(mu[oth]) * phi,
                            (1 - inv_logit(mu[oth])) * phi);
    return ll_zer + ll_one + ll_oth;
  }
