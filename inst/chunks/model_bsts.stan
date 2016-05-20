  loclev[1] ~ normal(Y[1], sigmaLL);
  for (n in 2:N) {
    if (tg[n] == tg[n - 1]) {
      loclev[n] ~ normal(loclev[n - 1], sigmaLL); 
    } else {
      loclev[n] ~ normal(Y[n], sigmaLL);
    }
  }
