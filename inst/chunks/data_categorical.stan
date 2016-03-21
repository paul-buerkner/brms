  int<lower=1> N_trait;
  int Y[N_trait];  // response variable
  int<lower=2> ncat;  // number of categories
  int J_trait[N_trait, ncat - 1];  // helps with indexing eta;
  // int J_trait[N_trait, 2];  // helps with indexing eta;
