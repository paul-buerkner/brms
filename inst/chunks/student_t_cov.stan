  /* multi-student-t log-PDF for special residual covariance structures 
   * currently only ARMA effects of order 1 are implemented 
   * Args: 
   *   y: response vector 
   *   nu: degrees of freedom parameter 
   *   eta: linear predictor 
   *   squared_se: square of the user defined standard errors 
   *               will be set to zero if non are defined 
   *   N_tg: number of groups 
   *   begin: indicates the first observation in each group 
   *   nrows: number of observations in each group 
   *   res_cov_matrix: AR1, MA1, or ARMA1 covariance matrix; 
   * Returns: 
   *   sum of the log-PDF values of all observations 
   */ 
   real student_t_cov_log(vector y, real nu, vector eta,  
                          vector squared_se, int N_tg, int[] begin, 
                          int[] nrows, matrix res_cov_matrix) { 
     vector[N_tg] log_post; 
     for (i in 1:N_tg) { 
       matrix[nrows[i], nrows[i]] Sigma; 
       vector[nrows[i]] y_part; 
       vector[nrows[i]] eta_part; 
       vector[nrows[i]] squared_se_part; 
       y_part <- segment(y, begin[i], nrows[i]); 
       eta_part <- segment(eta, begin[i], nrows[i]); 
       squared_se_part <- segment(squared_se, begin[i], nrows[i]); 
       Sigma <- block(res_cov_matrix, 1, 1, nrows[i], nrows[i]) 
                + diag_matrix(squared_se_part); 
       log_post[i] <- multi_student_t_log(y_part, nu, eta_part, Sigma); 
     }                        
     return sum(log_post); 
   }
