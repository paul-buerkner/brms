  /* dirichlet-multinomial-logit log-PDF
   * Args:
   *   y: vector of real response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real dirichlet_multinomial_logit2_lpmf(array[] int y, vector mu, real phi) {
      real alpha_plus = sum(phi * softmax(mu));
      return lgamma(alpha_plus) + sum(lgamma(phi * softmax(mu) + to_vector(y))) - lgamma(alpha_plus+sum(y)) - sum(lgamma(phi * softmax(mu)));
  }
