  // Explicit substitution avoids general matrix solves for one or two rows.
  vector mdivide_left_tri_low_brms(matrix L, vector b) {
    int K = rows(L);
    if (cols(L) != K || num_elements(b) != K) {
      return mdivide_left_tri_low(L, b);
    }
    if (K == 1 && L[1, 1] != 0.0) {
      return b / L[1, 1];
    }
    if (K == 2 && L[1, 1] != 0.0 && L[2, 2] != 0.0) {
      vector[2] x;
      x[1] = b[1] / L[1, 1];
      x[2] = (b[2] - L[2, 1] * x[1]) / L[2, 2];
      return x;
    }
    return mdivide_left_tri_low(L, b);
  }

  matrix mdivide_left_tri_low_brms(matrix L, matrix B) {
    int K = rows(L);
    if (cols(L) != K || rows(B) != K) {
      return mdivide_left_tri_low(L, B);
    }
    if (K == 1 && L[1, 1] != 0.0) {
      return B / L[1, 1];
    }
    if (K == 2 && L[1, 1] != 0.0 && L[2, 2] != 0.0) {
      matrix[2, cols(B)] X;
      X[1] = B[1] / L[1, 1];
      X[2] = (B[2] - L[2, 1] * X[1]) / L[2, 2];
      return X;
    }
    return mdivide_left_tri_low(L, B);
  }

  row_vector mdivide_right_tri_low_brms(row_vector b, matrix L) {
    int K = rows(L);
    if (cols(L) != K || num_elements(b) != K) {
      return mdivide_right_tri_low(b, L);
    }
    if (K == 1 && L[1, 1] != 0.0) {
      return b / L[1, 1];
    }
    if (K == 2 && L[1, 1] != 0.0 && L[2, 2] != 0.0) {
      row_vector[2] x;
      x[2] = b[2] / L[2, 2];
      x[1] = (b[1] - x[2] * L[2, 1]) / L[1, 1];
      return x;
    }
    return mdivide_right_tri_low(b, L);
  }

  matrix cholesky_decompose_brms(matrix A) {
    int K = rows(A);
    if (cols(A) != K) {
      return cholesky_decompose(A);
    }
    if (K == 1 && A[1, 1] > 0.0 && !is_inf(A[1, 1])) {
      return rep_matrix(sqrt(A[1, 1]), 1, 1);
    }
    if (K == 2 && A[1, 2] == A[2, 1] && A[1, 1] > 0.0 &&
        !is_inf(A[1, 1]) && !is_inf(A[2, 1]) && !is_inf(A[2, 2])) {
      matrix[2, 2] L = rep_matrix(0.0, 2, 2);
      real pivot;
      L[1, 1] = sqrt(A[1, 1]);
      L[2, 1] = A[2, 1] / L[1, 1];
      pivot = A[2, 2] - square(L[2, 1]);
      if (pivot > 0.0) {
        L[2, 2] = sqrt(pivot);
        return L;
      }
    }
    // Preserve Stan's validation and treatment of numerical boundary cases.
    return cholesky_decompose(A);
  }

  matrix chol2inv_brms(matrix L) {
    int K = rows(L);
    if (cols(L) != K) {
      return chol2inv(L);
    }
    if (K == 1 && L[1, 1] != 0.0) {
      return rep_matrix(inv_square(L[1, 1]), 1, 1);
    }
    if (K == 2 && L[1, 2] == 0.0 &&
        L[1, 1] != 0.0 && L[2, 2] != 0.0) {
      matrix[2, 2] precision;
      real a = inv(L[1, 1]);
      real d = inv(L[2, 2]);
      real c = -L[2, 1] * a / L[2, 2];
      precision[1, 1] = square(a) + square(c);
      precision[1, 2] = c * d;
      precision[2, 1] = precision[1, 2];
      precision[2, 2] = square(d);
      return precision;
    }
    return chol2inv(L);
  }
