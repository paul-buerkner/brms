  /* compute the kronecker product
   * Args: 
   *   A,B: matrices 
   * Returns: 
   *   kronecker product of A and B
   */ 
  matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:cols(A)) { 
      for (j in 1:rows(A)) { 
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
      } 
    } 
    return kron; 
  } 
  /* compute the kronecker product of two cholesky factors 
   * and multiply it by a vector, that is, (LA %kron% LB) * x
   * Args: 
   *   LA, LB: cholesky factors of matrices
   *   x: real vector
   * Returns: 
   *   a real vector of the same length as x
   */
  vector chol_kronecker_multiply(matrix LA, matrix LB, vector x) {
    vector[num_elements(x)] out = rep_vector(0, num_elements(x));
    int cols_LB = cols(LB);
    for (iA in 1:cols(LA)) {
      for (jA in 1:iA) {
        if (LA[iA, jA] > 1e-10) { 
          // avoid calculating products between unrelated individuals
          for (iB in 1:cols(LB)) {
            for (jB in 1:iB){
              int k = (cols_LB * (iA - 1)) + iB;
              int l = (cols_LB * (jA - 1)) + jB;
              out[k] = out[k] + LA[iA, jA] * LB[iB, jB] * x[l];
            }
          }
        }
      }
    }
    return out;
  }
