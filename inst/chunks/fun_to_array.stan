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
      Y[i] <- X[((i - 1) * K + 1):(i * K)]; 
    return Y; 
  } 
