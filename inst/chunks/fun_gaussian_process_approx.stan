  /* Spectral density function of a Gaussian process
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   * Returns: 
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp_quad(vector[] x, real sdgp, real lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    real constant = square(sdgp) * (sqrt(2 * pi()) * lscale)^D;
    real neg_half_lscale2 = -0.5 * square(lscale);
    vector[NB] out;
    for (m in 1:NB) {
      out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
    }
    return out;
  }
  /* compute an approximate 1D latent Gaussian process
   * Args:
   *   X: Matrix of Laplacian eigen functions at the covariate values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   *   zgp: vector of independent standard normal variables 
   *   slambda: square root of the Laplacian eigen values
   * Returns:  
   *   a vector to be added to the linear predictor
   */ 
  vector gpa(matrix X, real sdgp, real lscale, vector zgp, vector[] slambda) { 
    vector[cols(X)] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
    return X * (diag_spd .* zgp);
  }
