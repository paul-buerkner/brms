  // are two 1D integer arrays equal?
  int is_equal(array[] int a, array[] int b) {
    int n_a = size(a);
    int n_b = size(b);
    if (n_a != n_b) {
      return 0;
    }
    for (i in 1:n_a) {
      if (a[i] != b[i]) {
        return 0;
      }
    }
    return 1;
  }

