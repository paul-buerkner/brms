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
