  int<lower=1> N_trait;  // number of observations per category
  int Y[N_trait];  // response variable
  int<lower=2> ncat;  // number of categories
  int J_trait[N_trait, ncat - 1];  // helps with indexing eta;
