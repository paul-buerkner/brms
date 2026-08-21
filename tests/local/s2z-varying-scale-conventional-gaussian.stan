data {
  int<lower=1> N;
  int<lower=2> J;
  int<lower=2> M;
  vector[N] y;
  matrix[N, M - 1] Xc;
  matrix[N, M] Z;
  vector[M - 1] means_X;
  array[N] int<lower=1, upper=J> group;
  real<lower=0> beta_prior_sd;
  vector<lower=0>[M] sd_prior_rate;
  vector<lower=0>[M] sdlog_prior_sd;
  real<lower=0> sigma_prior_rate;
  real<lower=0> lkj_eta;
}
parameters {
  // brms places the intercept prior on the intercept for centered X.
  vector[M] beta_centered;
  real<lower=0> sigma;
  vector<lower=0>[M] sd;
  vector<lower=0>[M] sdlog;
  cholesky_factor_corr[M] L;
  matrix[M, J] z_effect;
  matrix[M, J] z_scale;
}
transformed parameters {
  matrix<lower=0>[J, M] sd_level;
  matrix[J, M] r;
  vector[N] mu;

  for (j in 1 : J) {
    vector[M] scale_j;
    matrix[M, M] L_j;
    for (k in 1 : M) {
      scale_j[k] = sd[k] * exp(sdlog[k] * z_scale[k, j]);
      sd_level[j, k] = scale_j[k];
    }
    L_j = diag_pre_multiply(scale_j, L);
    r[j] = (L_j * z_effect[, j])';
  }

  mu = rep_vector(beta_centered[1], N) +
    Xc * tail(beta_centered, M - 1);
  for (n in 1 : N) {
    mu[n] += dot_product(Z[n], r[group[n]]);
  }
}
model {
  beta_centered ~ normal(0, beta_prior_sd);
  sigma ~ exponential(sigma_prior_rate);
  sd ~ exponential(sd_prior_rate);
  sdlog ~ normal(0, sdlog_prior_sd);
  L ~ lkj_corr_cholesky(lkj_eta);
  to_vector(z_effect) ~ std_normal();
  to_vector(z_scale) ~ std_normal();
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[M] b_public;
  corr_matrix[M] Cor = multiply_lower_tri_self_transpose(L);
  matrix[J, M] finite_coefficient;

  b_public[1] = beta_centered[1] -
    dot_product(means_X, tail(beta_centered, M - 1));
  b_public[2 : M] = beta_centered[2 : M];
  for (j in 1 : J) {
    finite_coefficient[j] = b_public' + r[j];
  }
}
