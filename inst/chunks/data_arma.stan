  // data needed for ARMA effects 
  int<lower=0> Kar;  // AR order 
  int<lower=0> Kma;  // MA order 
  int<lower=1> Karma;  // max(Kma, Kar) 
  vector[N] tg;  // indicates independent groups
