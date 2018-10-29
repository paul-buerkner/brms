  /* Spectral density function of a 1D Gaussian process
   * Args:
   *   x: vector of numeric values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   * Returns: 
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp_quad(vector x, real sdgp, real lscale) {
    vector[rows(x)] out;
    for (m in 1:rows(x)) {
      out[m] = sdgp^2 * sqrt(2 * pi()) * 
        lscale * exp(-0.5 * lscale^2 * x[m]^2);
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
  vector gpa(matrix X, real sdgp, real lscale, vector zgp, vector slambda) { 
    vector[cols(X)] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
    return X * (diag_spd .* zgp);
  }
