  // Mode parameterization: w(k) = exp(nu * (k * log_mu - lgamma(k + 1))).
  // Sum relative to the modal term in both directions. The geometrically
  // bounded omitted mass AND first moment must meet the relative tolerance.
  // M limits work per direction, not the support of the distribution.
  vector com_poisson_sum(real log_mu, real nu, int M, real tol) {
    real mode = floor(exp(log_mu));
    real mass = 1;
    real log_moment = mode > 0 ? log(mode) : negative_infinity();
    vector[2] out;
    if (nu <= 0 || is_inf(nu) || is_nan(nu) || is_inf(log_mu) ||
        is_nan(log_mu) || is_inf(mode)) {
      reject("COM-Poisson requires finite positive mu and shape");
    }
    if (M < 0 || tol <= 0 || tol >= 1) reject("Invalid COM-Poisson summation controls");
    for (side in 1:2) {
      real k = mode;
      real lt = 0;
      int direction = side == 1 ? -1 : 1;
      int steps = 0;
      while (1) {
        real lr;
        real tail = positive_infinity();
        real moment_tail = positive_infinity();
        if (direction == -1 && k == 0) break;
        lr = direction == -1 ? nu * (log(k) - log_mu) :
                               nu * (log_mu - log(k + 1));
        // Ratios decrease away from the mode, so the next ratio bounds the
        // remainder. Split the error budget equally between the two tails.
        if (lr < 0) {
          real denom = -expm1(lr);
          real factor = direction == -1 ? fmax(k - 1, 0) :
                                          k + 1 + exp(lr) / denom;
          tail = lt + lr - log(denom);
          moment_tail = factor > 0 ? tail + log(factor) : negative_infinity();
        }
        if (tail <= log(tol / 2) + log(mass) &&
            moment_tail <= log(tol / 2) + log_moment) break;
        if (steps == M) reject("COM-Poisson summation failed to converge within M terms per direction");
        lt += lr;
        k += direction;
        mass += exp(lt);
        if (k > 0) log_moment = log_sum_exp(log_moment, log(k) + lt);
        steps += 1;
      }
    }
    out[1] = nu * (mode * log_mu - lgamma(mode + 1)) + log(mass);
    out[2] = exp(log_moment - log(mass));
    return out;
  }

  real log_Z_com_poisson(real log_mu, real nu) {
    return com_poisson_sum(log_mu, nu, 10000, 1e-12)[1];
  }

  // Do not special-case nu == 1: the value is Poisson there but the derivative
  // with respect to a free shape parameter is not zero.
  real com_poisson_log_lpmf(int y, real log_mu, real nu) {
    if (y < 0) return negative_infinity();
    return nu * (y * log_mu - lgamma(y + 1.0)) - log_Z_com_poisson(log_mu, nu);
  }
  real com_poisson_lpmf(int y, real mu, real nu) {
    return com_poisson_log_lpmf(y | log(mu), nu);
  }

  // Unnormalized one-sided sum. Log arithmetic also allows crossing the mode
  // if the first tail selected by the CDF turns out to have probability > 1/2.
  real com_poisson_log_tail(int y, real log_mu, real nu, int lower_tail) {
    real k = lower_tail ? y : y + 1.0;
    real first = nu * (k * log_mu - lgamma(k + 1));
    real total = 0;
    real lt = 0;
    int steps = 0;
    while (1) {
      real lr;
      if (lower_tail && k == 0) break;
      lr = lower_tail ? nu * (log(k) - log_mu) : nu * (log_mu - log(k + 1));
      if (lr < 0 && lt + lr - log(-expm1(lr)) <= log(1e-12) + total) break;
      if (steps == 10000) reject("COM-Poisson tail sum failed to converge within M terms");
      lt += lr;
      k += lower_tail ? -1 : 1;
      total = log_sum_exp(total, lt);
      steps += 1;
    }
    return first + total;
  }

  real com_poisson_log_cumulative(int y, real mu, real nu, int lower_tail) {
    real log_mu;
    real log_Z;
    real lp;
    int direct_lower;
    if (y < 0) return lower_tail ? negative_infinity() : 0;
    log_mu = log(mu);
    log_Z = log_Z_com_poisson(log_mu, nu);
    direct_lower = y < floor(mu);
    lp = com_poisson_log_tail(y, log_mu, nu, direct_lower) - log_Z;
    if (lp > -log(2)) {
      direct_lower = !direct_lower;
      lp = com_poisson_log_tail(y, log_mu, nu, direct_lower) - log_Z;
    }
    lp = fmin(lp, 0);
    if (direct_lower == lower_tail) return lp;
    if (lp == 0) return negative_infinity();
    return log1m_exp(lp);
  }
  real com_poisson_lcdf(int y, real mu, real nu) {
    return com_poisson_log_cumulative(y, mu, nu, 1);
  }
  real com_poisson_lccdf(int y, real mu, real nu) {
    return com_poisson_log_cumulative(y, mu, nu, 0);
  }
