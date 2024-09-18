  /* Spectral density function of a Gaussian process with Matern 3/2 kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric vector of length NB of the SPD evaluated at 'x'
   */
  vector spd_gp_matern32(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    real constant = square(sdgp) *
      (2^D * pi()^(D / 2.0) * tgamma((D + 3.0) / 2) * 3^(3.0 / 2)) /
      (0.5 * sqrt(pi()));
    real expo = -(D + 3.0) / 2;
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real lscale2 = square(lscale[1]);
      constant = constant * lscale[1]^D;
      for (m in 1:NB) {
        out[m] = constant * (3 + lscale2 * dot_self(x[m]))^expo;
      }
    } else {
      // multi-dimensional non-isotropic GP
      vector[Dls] lscale2 = square(lscale);
      constant = constant * prod(lscale);
      for (m in 1:NB) {
        out[m] = constant * (3 + dot_product(lscale2, square(x[m])))^expo;
      }
    }
    return out;
  }
