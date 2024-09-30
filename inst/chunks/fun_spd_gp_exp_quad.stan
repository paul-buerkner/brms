  /* Spectral density function of a Gaussian process
   * with squared exponential covariance kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric vector of length NB of the SPD evaluated at 'x'
   */
  vector spd_gp_exp_quad(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    real constant = square(sdgp) * sqrt(2 * pi())^D;
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real neg_half_lscale2 = -0.5 * square(lscale[1]);
      constant = constant * lscale[1]^D;
      for (m in 1:NB) {
        out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
      }
    } else {
      // multi-dimensional non-isotropic GP
      vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
      constant = constant * prod(lscale);
      for (m in 1:NB) {
        out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
      }
    }
    return out;
  }
