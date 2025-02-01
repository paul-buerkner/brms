  /* dirichlet-multinomial-logit log-PMF
   * Args:
   *   y: array of integer response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter (sum of Dirichlet alphas)
   * Returns:
   *   a scalar to be added to the log posterior
  */
  real dirichlet_multinomial_logit2_lpmf(array[] int y, vector mu, real phi) {
    // get Dirichlet alphas
    int N = num_elements(mu);
    vector[N] alpha = phi * softmax(mu);

    // get trials from y
    real T = sum(y);

    return lgamma(phi) + lgamma(T + 1.0) - lgamma(T + phi) +
      sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha)) - sum(lgamma(to_vector(y) + 1));
  }
