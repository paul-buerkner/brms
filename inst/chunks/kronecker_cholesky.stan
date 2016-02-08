  /* calculate the cholesky factor of a kronecker covariance matrix 
   * Args: 
   *   X: a covariance matrix 
   *   L: cholesky factor of another covariance matrix 
   *   sd: standard deviations for scaling 
   * Returns: 
   *   cholesky factor of kronecker(X, L * L')  
   */ 
  matrix kronecker_cholesky(matrix X, matrix L, vector sd) { 
    matrix[rows(X)*rows(L), cols(X)*cols(L)] kron; 
    matrix[rows(L), cols(L)] C; 
    int rX; 
    int rC; 
    C <- multiply_lower_tri_self_transpose(L); 
    rX <- rows(X); 
    rC <- rows(C); 
    for (i in 1:rX) { 
      for (j in 1:rC) { 
        for (k in 1:rX) { 
          for (l in 1:rC) { 
            kron[(k-1) * rC+l, (i-1) * rC+j] <- sd[l] * sd[j] * X[k,i] * C[l,j]; 
          } 
        } 
      } 
    } 
    return cholesky_decompose(kron); 
  } 
  /* turn a vector into a 2 dimensional array 
   * Args: 
   *   X: a vector 
   *   N: first dimension of the desired array 
   *   K: second dimension of the desired array 
   * Returns: 
   *   an array of dimension N x K 
   */ 
  vector[] to_array(vector X, int N, int K) { 
    vector[K] Y[N]; 
    for (i in 1:N) 
      Y[i] <- segment(X, (i - 1) * K + 1, K); 
    return Y; 
  } 
