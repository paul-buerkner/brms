  // The previous brms normalizer was inspired by code of Ben Goodrich and
  // improved following suggestions of Sebastian Weber (#892).
  // Mode-centred sum with geometric bounds on omitted mass; see Shmueli et al.
  // (2005), Appendix B.1, Eqs. (33)-(35), doi:10.1111/j.1467-9876.2005.00474.x.
  // Apply the bounds in both directions and additionally bound the first moment.
  // M limits terms per direction; tol controls relative mass/moment tails.
  // If the first omitted weight is a at count t and outward ratios are <= r < 1:
  //   omitted mass <= a / (1 - r)
  //   right first moment <= a * (t / (1 - r) + r / (1 - r)^2)
  //   left first moment <= t * a / (1 - r), since remaining counts are <= t.
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

  // Same empirical region as R: mu * nu >= 200 and mu / nu >= 50.
  // This is a tested approximation rule, not a uniform remainder bound.
  int com_poisson_use_approx(real log_mu, real nu) {
    return log_mu + log(nu) >= log(200) && log_mu - log(nu) >= log(50);
  }

  // Gaunt et al. (2019), Eq. (A.31), doi:10.1007/s10463-017-0629-6.
  // Retain the existing three corrections and derive the mean from them.
  vector com_poisson_approx(real log_mu, real nu) {
    real mu = exp(log_mu);
    real a = nu / mu;
    real b = 1 / (nu * mu);
    real t1 = (a - b) / 24;
    real t2 = (a - b) * (a + 23 * b) / 1152;
    real t3 = (a - b) * (5 * square(a) - 298 * a * b + 11237 * square(b)) / 414720;
    real correction = t1 + t2 + t3;
    vector[2] out;
    out[1] = nu * mu - (nu - 1) / 2 * (log(2 * pi()) + log_mu) -
             log(nu) / 2 + log1p(correction);
    out[2] = mu - (1 - 1 / nu) / 2 -
             (t1 + 2 * t2 + 3 * t3) / (nu * (1 + correction));
    return out;
  }

  vector com_poisson_normalizer(real log_mu, real nu) {
    if (nu <= 0 || is_inf(nu) || is_nan(nu) || is_inf(log_mu) ||
        is_nan(log_mu) || is_inf(exp(log_mu))) {
      reject("COM-Poisson requires finite positive mu and shape");
    }
    if (com_poisson_use_approx(log_mu, nu)) return com_poisson_approx(log_mu, nu);
    return com_poisson_sum(log_mu, nu, 10000, 1e-12);
  }

  real log_Z_com_poisson(real log_mu, real nu) {
    return com_poisson_normalizer(log_mu, nu)[1];
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
