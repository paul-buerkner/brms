  // log approximate normalizing constant of the COM poisson distribuion
  // based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
  // Args:
  //   log_lambda: log of the classical rate parameter, that is nu * log(mu)
  //   nu: positive shape parameter
  real log_Z_com_poisson_approx(real log_lambda, real nu) {
    real nu2 = nu^2;
    real log_common = log(nu) + log_lambda / nu;
    array[4] real resids;
    real ans;
    real lcte = (nu * exp(log_lambda / nu)) -
      ((nu - 1) / (2 * nu) * log_lambda +
       (nu - 1) / 2 * log(2 * pi()) + 0.5 * log(nu));
    real c_1 = (nu2 - 1) / 24;
    real c_2 = (nu2 - 1) / 1152 * (nu2 + 23);
    real c_3 = (nu2 - 1) / 414720 * (5 * square(nu2) - 298 * nu2 + 11237);
    resids[1] = 1;
    resids[2] = c_1 * exp(-1 * log_common);
    resids[3] = c_2 * exp(-2 * log_common);
    resids[4] = c_3 * exp(-3 * log_common);
    ans = lcte + log(sum(resids));
    return ans;
  }

  // log of kth term of the normalizing series of the COM Poisson distribution
  // Args:
  //   log_lambda: log of the classical rate parameter, that is nu * log(mu)
  //   nu: positive shape parameter
  //   k: k-th term
  real log_k_term(real log_lambda, real nu, int k) {
    return (k - 1) * log_lambda - nu * lgamma(k);
  }

  // bound for the remainder of the normalizing series of the COM Poisson
  // distribution given the last two terms in log-scale
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  real bound_remainder(real k_current_term, real k_previous_term) {
    return k_current_term - log(- expm1(k_current_term - k_previous_term));
  }

  // stopping criterio with bucket
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  //   k: k term of the series
  //   leps: log(eps)
  int stopping_criterio_bucket(real k_current_term, real k_previous_term, int k, real leps) {
    if (k % 50 == 0) {
      return (bound_remainder(k_current_term, k_previous_term) >= leps);
    }
    return (1e300 >= leps); // Int > leps
  }

  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // improved following suggestions of Sebastian Weber (#892)
  // brms parameterizes com_poisson through the modal location mu, while the
  // series below is evaluated in the classical parameterization lambda = mu^nu;
  // convert internally to log_lambda = nu * log_mu
  // Args:
  //   log_mu: log of the modal location parameter
  //   nu: positive shape parameter
  real log_Z_com_poisson(real log_mu, real nu) {
    real log_lambda = nu * log_mu;
    real log_Z;
    int k = 2;
    int M = 10000;
    real leps = -23 * log2();
    vector[M] log_Z_terms;

    if (nu == 1) {
      return exp(log_mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }
    // use the same rule-of-thumb region as log_Z_com_poisson() in R:
    // lambda >= 1.5 and mu = lambda^(1 / nu) >= 1.5
    if (log_lambda >= log(1.5) && log_mu >= log(1.5)) {
      return log_Z_com_poisson_approx(log_lambda, nu);
    }
    // direct computation of the truncated series
    // check if the Mth term of the series pass in the stopping criteria
    if (bound_remainder(log_k_term(log_lambda, nu, M),
                        log_k_term(log_lambda, nu, M - 1)) >= leps) {
      reject("nu is too close to zero.");
    }

    // first 2 terms of the series
    log_Z_terms[1] = log_k_term(log_lambda, nu, 1);
    log_Z_terms[2] = log_k_term(log_lambda, nu, 2);

    while (((log_Z_terms[k] >= log_Z_terms[k-1]) ||
      (stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps))) &&
      k < M) {
      k += 1;
      log_Z_terms[k] = log_k_term(log_lambda, nu, k);
    }
    log_Z = log_sum_exp(log_Z_terms[1:k]);

    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args:
  //   y: the response value
  //   log_mu: log of the modal location parameter
  //   nu: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_mu, real nu) {
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return nu * (y * log_mu - lgamma(y + 1))
           - log_Z_com_poisson(log_mu, nu);
  }
  // COM Poisson log-PMF for a single response
  real com_poisson_lpmf(int y, real mu, real nu) {
    if (nu == 1) return poisson_lpmf(y | mu);
    return com_poisson_log_lpmf(y | log(mu), nu);
  }
  // COM Poisson log-CDF for a single response
  real com_poisson_lcdf(int y, real mu, real nu) {
    real log_mu = log(mu);
    real log_lambda = nu * log_mu;
    real log_Z;  // log denominator
    vector[y] log_num_terms;  // terms of the log numerator
    if (nu == 1) {
      return poisson_lcdf(y | mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }
    if (y > 10000) {
      reject("cannot handle y > 10000");
    }
    if (y * log_lambda - nu * lgamma(y + 1) <= -36.0) {
      // y is large enough for the CDF to be very close to 1;
      return 0;
    }
    log_Z = log_Z_com_poisson(log_mu, nu);
    if (y == 0) {
      return -log_Z;
    }
    // first 2 terms of the series
    log_num_terms[1] = log1p_exp(log_lambda);
    // remaining terms of the series until y
    for (k in 2:y) {
      log_num_terms[k] = k * log_lambda - nu * lgamma(k + 1);
    }
    return log_sum_exp(log_num_terms) - log_Z;
  }
  // COM Poisson log-CCDF for a single response
  real com_poisson_lccdf(int y, real mu, real nu) {
    return log1m_exp(com_poisson_lcdf(y | mu, nu));
  }
