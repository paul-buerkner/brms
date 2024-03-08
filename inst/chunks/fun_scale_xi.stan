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
    vector[rows(y)] x = (y - mu) / sigma;
    vector[2] bounds = [-inv(min(x)), -inv(max(x))]';
    real lb = min(bounds);
    real ub = max(bounds);
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
  real scale_xi(real xi, vector y, vector mu, vector sigma) {
    vector[rows(y)] x = (y - mu) ./ sigma;
    vector[2] bounds = [-inv(min(x)), -inv(max(x))]';
    real lb = min(bounds);
    real ub = max(bounds);
    return inv_logit(xi) * (ub - lb) + lb;
  }
