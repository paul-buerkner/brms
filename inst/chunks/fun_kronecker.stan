  /* compute the kronecker product
   * Args: 
   *   A,B: matrices 
   * Returns: 
   *   kronecker product of A and B
   */ 
  matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:rows(A)) { 
      for (j in 1:rows(B)) { 
        for (k in 1:rows(A)) { 
          for (l in 1:rows(B)) { 
            kron[(k-1) * rows(B)+l, (i-1) * rows(B)+j] <- A[k,i] * B[l,j]; 
          } 
        } 
      } 
    } 
    return kron; 
  } 
