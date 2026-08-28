  vector sum_to_zero_constrain_brms(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1:N) {
      int i = N - ii + 1;
      real w = y[i] * inv_sqrt(i * (i + 1.0));
      sum_w += w;
      z[i] += sum_w;
      z[i + 1] -= i * w;
    }
    return z;
  }

  real s2z_require_finite_brms(real x) {
    if (is_nan(x) || is_inf(x)) {
      reject("S2Z population-prior locations must be finite.");
    }
    return x;
  }

  real s2z_require_positive_brms(real x) {
    if (is_nan(x) || is_inf(x) || x <= 0) {
      reject("S2Z population-prior scales and degrees of freedom must be ",
             "finite and strictly positive.");
    }
    return x;
  }

  real s2z_prior_coordinate_brms(real x, int index, int expected_size) {
    return x;
  }

  real s2z_prior_coordinate_brms(vector x, int index, int expected_size) {
    if (num_elements(x) != expected_size) {
      reject("An S2Z vector-valued population-prior argument must have one ",
             "entry per population-level coefficient.");
    }
    return x[index];
  }

  real s2z_prior_coordinate_brms(row_vector x, int index,
                                 int expected_size) {
    if (num_elements(x) != expected_size) {
      reject("An S2Z vector-valued population-prior argument must have one ",
             "entry per population-level coefficient.");
    }
    return x[index];
  }

  real s2z_prior_coordinate_brms(array[] real x, int index,
                                 int expected_size) {
    if (num_elements(x) != expected_size) {
      reject("An S2Z vector-valued population-prior argument must have one ",
             "entry per population-level coefficient.");
    }
    return x[index];
  }

  real s2z_prior_coordinate_brms(array[] int x, int index,
                                 int expected_size) {
    if (num_elements(x) != expected_size) {
      reject("An S2Z vector-valued population-prior argument must have one ",
             "entry per population-level coefficient.");
    }
    return x[index];
  }
