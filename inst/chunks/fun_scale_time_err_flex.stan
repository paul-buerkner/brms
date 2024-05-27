  /* scale and correlate time-series residuals
   * allows for flexible correlation matrix subsets
   * Deviating Args:
   *   Jtime: array of time indices per group
   * Returns:
   *   vector of scaled and correlated residuals
   */
   vector scale_time_err_flex(vector zerr, real sderr, matrix chol_cor,
                              array[] int nobs, array[] int begin,
                              array[] int end, array[,] int Jtime) {
     vector[rows(zerr)] err;
     int I = size(nobs);
     array[I] int has_err = rep_array(0, I);
     int i = 1;
     matrix[rows(chol_cor), cols(chol_cor)] L;
     matrix[rows(chol_cor), cols(chol_cor)] Cov;
     L = sderr * chol_cor;
     Cov = multiply_lower_tri_self_transpose(L);
     while (i <= I) {
       array[nobs[i]] int iobs = Jtime[i, 1:nobs[i]];
       matrix[nobs[i], nobs[i]] L_i;
       if (is_equal(iobs, sequence(1, rows(L)))) {
         // all timepoints are present in this group
         L_i = L;
       } else {
         // arbitrary subsets cannot be taken on L directly
         L_i = cholesky_decompose(Cov[iobs, iobs]);
       }
       err[begin[i]:end[i]] = L_i * zerr[begin[i]:end[i]];
       has_err[i] = 1;
       // find all additional groups where we have the same timepoints
       for (j in (i+1):I) {
         if (has_err[j] == 0 && is_equal(Jtime[j], Jtime[i]) == 1) {
           err[begin[j]:end[j]] = L_i * zerr[begin[j]:end[j]];
           has_err[j] = 1;
         }
       }
       while (i <= I && has_err[i] == 1) {
         i += 1;
       }
    }
    return err;
  }
