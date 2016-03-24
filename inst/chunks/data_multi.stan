  int<lower=1> N_trait;  // number of observations per trait 
  int<lower=1> K_trait;  // number of responses   
  int NC_trait;  // number of residual correlations 
  vector[K_trait] Y[N_trait];  // response matrix
