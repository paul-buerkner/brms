  /* scale auxiliary parameter xi to a suitable region
   * expecting sigma to be a scalar
   * Args:
   *   xi: unscaled shape parameter
   *   y: response values
   *   mu: location parameter
   *   sigma: scale parameter
   * Returns:
   *   scaled shape parameter xi
   */
  real scale_xi(real xi, vector y, vector mu, real sigma) {
    vector[2] bounds;
    real lb;
    real ub;
    bounds[1] = - inv(min((y - mu) / sigma));
    bounds[2] = - inv(max((y - mu) / sigma));
    lb = min(bounds);
    ub = max(bounds);
    return inv_logit(xi) * (ub - lb) + lb;
  }
  /* scale auxiliary parameter xi to a suitable region
   * expecting sigma to be a vector
   * Args:
   *   xi: unscaled shape parameter
   *   y: response values
   *   mu: location parameter
   *   sigma: scale parameter
   * Returns:
   *   scaled shape parameter xi
   */
  real scale_xi_vector(real xi, vector y, vector mu, vector sigma) {
    vector[2] bounds;
    real lb;
    real ub;
    bounds[1] = - inv(min((y - mu) ./ sigma));
    bounds[2] = - inv(max((y - mu) ./ sigma));
    lb = min(bounds);
    ub = max(bounds);
    return inv_logit(xi) * (ub - lb) + lb;
  }
