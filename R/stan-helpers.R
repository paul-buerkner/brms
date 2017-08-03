stan_autocor <- function(autocor, bterms, family, prior) {
  # Stan code related to autocorrelation structures
  # Args:
  #   autocor: object of class cor_brms
  #   bterms: object of class brmsterms
  #   family: the model family
  #   prior: object of class brmsprior
  stopifnot(is.family(family))
  stopifnot(is.brmsterms(bterms))
  is_linear <- is_linear(family)
  resp <- bterms$response
  is_mv <- is_linear && length(resp) > 1L
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  Karr <- get_arr(autocor)
  out <- list()
  if (Kar || Kma) {
    if (!is_linear) {
      stop2("ARMA models are not yet implemented ", 
            "for family '", family$family, "'.") 
    }
    str_add(out$data) <- paste0( 
      "  // data needed for ARMA correlations \n",
      "  int<lower=0> Kar;  // AR order \n",
      "  int<lower=0> Kma;  // MA order \n"
    )
    str_add(out$tdataD) <- paste0( 
      "  int max_lag = max(Kar, Kma); \n"
    )
    # restrict ARMA effects to be in [-1,1] when using covariance
    # formulation as they cannot be outside this interval anyway
    if (Kar) {
      ar_bound <- with(prior, bound[class == "ar"])
      str_add(out$par) <- paste0( 
        "  vector", ar_bound, "[Kar] ar;  // autoregressive effects \n"
      )
      str_add(out$prior) <- paste0(stan_prior(prior, class = "ar"))
    }
    if (Kma) {
      ma_bound <- with(prior, bound[class == "ma"])
      str_add(out$par) <- paste0( 
        "  vector", ma_bound, "[Kma] ma;  // moving-average effects \n"
      )
      str_add(out$prior) <- paste0(stan_prior(prior, class = "ma"))
    }
    if (use_cov(autocor)) {
      # if the user wants ARMA effects to be estimated using
      # a covariance matrix for residuals
      err_msg <- "ARMA covariance matrices are not yet working"
      if (is_mv) {
        stop2(err_msg, " in multivariate models.")
      }
      if (is.formula(bterms$adforms$disp)) {
        stop2(err_msg, " when specifying 'disp'.")
      }
      if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
        stop2(err_msg, " when predicting 'sigma' or 'nu'.")
      }
      str_add(out$data) <- paste0( 
        "  #include 'data_arma_cov.stan' \n"
      )
      if (!is.formula(bterms$adforms$se)) {
        str_add(out$tdataD) <- "  vector[N] se2 = rep_vector(0, N); \n"
      }
      str_add(out$tparD) <- paste0(
        "  matrix[max(nobs_tg), max(nobs_tg)] res_cov_matrix; \n"                  
      )
      if (Kar && !Kma) {
        cov_mat_fun <- "ar1"
        cov_mat_args <- "ar[1]"
      } else if (!Kar && Kma) {
        cov_mat_fun <- "ma1"
        cov_mat_args <- "ma[1]"
      } else {
        cov_mat_fun <- "arma1"
        cov_mat_args <- "ar[1], ma[1]"
      }
      str_add(out$tparC1) <- paste0(
        "  // compute residual covariance matrix \n",
        "  res_cov_matrix = cov_matrix_", cov_mat_fun, 
        "(", cov_mat_args, ", sigma, max(nobs_tg)); \n"
      )
      # defined self-made functions for the functions block
      if (family$family == "gaussian") {
        str_add(out$fun) <- paste0(
          "  #include 'fun_normal_cov.stan' \n"
        )
      } else {  
        str_add(out$fun) <- paste0(
          "  #include 'fun_student_t_cov.stan' \n"
        )
      }
      if (Kar && !Kma) {
        str_add(out$fun) <- paste0(
          "  #include 'fun_cov_matrix_ar1.stan' \n"
        )
      } else if (!Kar && Kma) {
        str_add(out$fun) <- paste0(
          "  #include 'fun_cov_matrix_ma1.stan' \n"
        )
      } else {
        str_add(out$fun) <- paste0(
          "  #include 'fun_cov_matrix_arma1.stan' \n"
        )
      }
    } else {
      err_msg <- "Please set cov = TRUE in cor_arma / cor_ar / cor_ma"
      if (is.formula(bterms$adforms$se)) {
        stop2(err_msg, " when specifying 'se'.")
      }
      if (length(bterms$dpars[["mu"]]$nlpars)) {
        stop2(err_msg, " in non-linear models.")
      }
      if (!identical(family$link, "identity")) {
        stop2(err_msg, " when using non-identity links.")
      }
      str_add(out$data) <- paste0( 
        "  int<lower=0> J_lag[N]; \n"                
      )
      if (is_mv) {
        rs <- usc(resp, "prefix")
        index <- paste0("n, ", seq_along(rs))
      } else {
        rs <- ""
        index <- "n"
      }
      str_add(out$modelD) <- paste0(
        "  // objects storing residuals \n",
        collapse(
          "  matrix[N, max_lag] E", rs, 
          " = rep_matrix(0, N, max_lag); \n",
          "  vector[N] e", rs, "; \n"
        )
      )
      str_add(out$modelC2) <- paste0(
        "    // computation of ARMA correlations \n",
        collapse(
          "    e", rs, "[n] = Y[", index, "] - mu", rs, "[n]", "; \n"
        ),
        "    for (i in 1:J_lag[n]) { \n",
        collapse(
         "      E", rs, "[n + 1, i] = e", rs, "[n + 1 - i]; \n"
        ),
        "    } \n"
      )
    } 
  }
  if (Karr) {
    # autoregressive effects of the response
    err_msg <- "ARR models are not yet working"
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " in non-linear models.")
    }
    if (is_mv) {
      stop2(err_msg, " in multivariate models.")
    }
    str_add(out$data) <- paste0(
      "  // data needed for ARR correlations \n",
      "  int<lower=1> Karr; \n",
      "  matrix[N, Karr] Yarr;  // ARR design matrix \n"
    )
    str_add(out$par) <- paste0(
      "  vector", with(prior, bound[class == "arr"]), "[Karr] arr;",
      "  // autoregressive effects of the response \n"
    )
    str_add(out$prior) <- paste0(stan_prior(prior, class = "arr"))
  }
  if (is.cor_sar(autocor)) {
    err_msg <- "SAR models are not yet working"
    if (!is_linear) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (is_mv) {
      stop2(err_msg, " in multivariate models.")
    }
    if (is.formula(bterms$adforms$disp)) {
      stop2(err_msg, " when specifying 'disp'.")
    }
    if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
      stop2(err_msg, " when predicting 'sigma' or 'nu'.")
    }
    str_add(out$data) <- paste0(
      "  matrix[N, N] W;  // spatial weight matrix \n"                  
    )
    if (identical(autocor$type, "lag")) {
      if (family$family == "gaussian") {
        str_add(out$fun) <- paste0(
          "  #include 'fun_normal_lagsar.stan' \n"
        ) 
      } else if (family$family == "student") {
        str_add(out$fun) <- paste0(
          "  #include 'fun_student_t_lagsar.stan' \n"
        )
      }
      str_add(out$par) <- paste0( 
        "  real<lower=0,upper=1>  lagsar;  // SAR parameter \n"
      )
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "lagsar")
      )
    } else if (identical(autocor$type, "error")) {
      if (family$family == "gaussian") {
        str_add(out$fun) <- paste0(
          "  #include 'fun_normal_errorsar.stan' \n"
        ) 
      } else if (family$family == "student") {
        str_add(out$fun) <- paste0(
          "  #include 'fun_student_t_errorsar.stan' \n"
        ) 
      }
      str_add(out$par) <- paste0( 
        "  real<lower=0,upper=1> errorsar;  // SAR parameter \n"
      )
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "errorsar")
      )
    }
  }
  if (is.cor_car(autocor)) {
    err_msg <- "CAR models are not yet working"
    if (is_mv) {
      stop2(err_msg, " in multivariate models.")
    }
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " in non-linear models.")
    }
    str_add(out$data) <- paste0(
      "  // data for the CAR structure \n",
      "  int<lower=1> Nloc; \n",
      "  vector[Nloc] Nneigh; \n",
      "  vector[Nloc] eigenW; \n",
      "  int<lower=1> Jloc[N]; \n",
      "  int<lower=0> Nedges; \n",
      "  int<lower=1> edges1[Nedges]; \n",
      "  int<lower=1> edges2[Nedges]; \n"
    )
    str_add(out$par) <- paste0(
      "  // parameters for the CAR structure \n",
      "  real<lower=0> sdcar; \n"
    )
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "sdcar")
    )
    if (identical(autocor$type, "escar")) {
      str_add(out$fun) <- paste0(
        "  #include 'fun_sparse_car_lpdf.stan' \n"        
      )
      str_add(out$par) <- paste0(
        "  real<lower=0, upper=1> car; \n",
        "  vector[Nloc] rcar; \n"
      )
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "car"),
        "  target += sparse_car_lpdf(\n", 
        "    rcar | car, sdcar, Nloc, Nedges,\n",
        "    Nneigh, eigenW, edges1, edges2\n", 
        "  ); \n"
      )
    } else if (identical(autocor$type, "esicar")) {
      str_add(out$fun) <- paste0(
        "  #include 'fun_sparse_icar_lpdf.stan' \n"        
      )
      str_add(out$par) <- paste0(
        "  vector[Nloc - 1] zcar; \n"
      )
      str_add(out$tparD) <- paste0(
        "  vector[Nloc] rcar; \n"                
      )
      str_add(out$tparC1) <- paste0(
        "  // apply sum-to-zero constraint \n",
        "  rcar[1:(Nloc - 1)] = zcar; \n",
        "  rcar[Nloc] = - sum(zcar); \n"
      )
      str_add(out$prior) <- paste0(
        "  target += sparse_icar_lpdf(\n", 
        "    rcar | sdcar, Nloc, Nedges,\n",
        "    Nneigh, eigenW, edges1, edges2\n", 
        "  ); \n"
      )
    } 
  }
  if (is.cor_bsts(autocor)) {
    err_msg <- "BSTS models are not yet working"
    if (is_ordinal(family) || 
        family$family %in% c("bernoulli", "categorical")) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (is_mv) {
      stop2(err_msg, " in multivariate models.")
    }
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " in non-linear models.")
    }
    str_add(out$data) <- 
      "  vector[N] tg;  // indicates independent groups \n"
    str_add(out$par) <- paste0(
      "  vector[N] loclev;  // local level terms \n",
      "  real<lower=0> sigmaLL;  // SD of local level terms \n"
    )
    if (is_linear && !is_mv) {
      # ensures that the means of the loclev priors start close 
      # to the response values; this often helps with convergence
      link <- stan_link(family$link)
      center <- paste0(link, "(Y[", c("1", "n"), "])")
    } else {
      center <- c("0", "0")
    }
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "sigmaLL"),
      "  target += normal_lpdf(loclev[1] | ", center[1], ", sigmaLL); \n",
      "  for (n in 2:N) { \n",
      "    if (tg[n] == tg[n - 1]) { \n",
      "      target += normal_lpdf(loclev[n] | loclev[n - 1], sigmaLL); \n",
      "    } else { \n",
      "      target += normal_lpdf(loclev[n] | ", center[2], ", sigmaLL); \n",
      "    } \n",
      "  } \n"
    )
  }
  if (is.cor_fixed(autocor)) {
    if (!is_linear) {
      stop2("Fixed residual covariance matrices are not yet ", 
            "implemented for family '", family$family, "'.") 
    }
    str_add(out$data) <- paste0( 
      "  matrix[N, N] V;  // known residual covariance matrix \n"
    )
    if (family$family %in% "gaussian") {
      str_add(out$tdataD) <- paste0(
        "  matrix[N, N] LV = cholesky_decompose(V); \n"
      )
    }
  }
  out
}

stan_mv <- function(family, response, prior) {
  # some Stan code for multivariate models
  # Args:
  #   family: model family
  #   response: names of the response variables
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # Returns: 
  #   list containing Stan code specific for multivariate models
  stopifnot(is.family(family))
  out <- list()
  nresp <- length(response)
  if (nresp > 1L) {
    if (is_linear(family)) {
      str_add(out$data) <- "  #include 'data_mv.stan' \n"
      str_add(out$par) <- paste0(
        "  // parameters for multivariate linear models \n",
        "  vector<lower=0>[nresp] sigma; \n",
        "  cholesky_factor_corr[nresp] Lrescor; \n"
      )
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "sigma", coef = response),
        stan_prior(prior, class = "Lrescor")
      )
      if (family$family == "gaussian") {
        str_add(out$tparD) <- paste0(
          "  // cholesky factor of residual covariance matrix \n",
          "  cholesky_factor_cov[nresp] LSigma = ",
          "diag_pre_multiply(sigma, Lrescor); \n"
        )
      } else if (family$family == "student") {
        str_add(out$tparD) <- paste0(
          "  // residual covariance matrix \n",
          "  cov_matrix[nresp] Sigma",
          " = multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor)); \n"
        )
      }
      str_add(out$genD) <- paste0(
        "  matrix[nresp, nresp] Rescor",
        " = multiply_lower_tri_self_transpose(Lrescor); \n",
        "  vector<lower=-1,upper=1>[nrescor] rescor; \n"
      )
      rescor_genC <- ulapply(2:nresp, function(i) 
        lapply(1:(i - 1), function(j) paste0(
          "  rescor[", (i - 1) * (i - 2) / 2 + j, 
          "] = Rescor[", j, ", ", i, "]; \n"
        ))
      )
      str_add(out$genC) <- paste0(
        "  // take only relevant parts of residual correlation matrix \n",
        collapse(rescor_genC)
      )
    } else if (!is_categorical(family)) {
      stop2("Multivariate models are not yet implemented ", 
            "for family '", family$family, "'.")
    }
  }
  out
}

stan_ordinal <- function(family, prior, cs, disc) {
  # Ordinal effects in Stan
  # Args:
  #   family: the model family
  #   prior: object of class brmsprior
  #   cs: logical; are there category specific effects?
  #   disc: logical; discrimination parameter used?
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  stopifnot(is.family(family))
  out <- list()
  if (is_ordinal(family)) {
    # define Stan code similar for all ordinal models
    str_add(out$data) <- "  int ncat;  // number of categories \n"
    th <- function(k, fam = family) {
      # helper function generating stan code inside ilink(.)
      sign <- ifelse(fam %in% c("cumulative", "sratio"), " - ", " + ")
      ptl <- ifelse(cs, paste0(sign, "mucs[k]"), "") 
      if (sign == " - ") {
        out <- paste0("thres[", k, "]", ptl, " - mu")
      } else {
        out <- paste0("mu", ptl, " - thres[", k, "]")
      }
      paste0("disc * (", out, ")")
    }
    link <- family$link
    threshold <- family$threshold
    family <- family$family
    ilink <- stan_ilink(link)
    type <- ifelse(family == "cumulative", "ordered", "vector")
    intercept <- paste0(
      "  ", type, "[ncat-1] temp_Intercept;",
      "  // temporary thresholds \n"
    )
    if (threshold == "flexible") {
      str_add(out$par) <- intercept
      str_add(out$prior) <- stan_prior(prior, class = "temp_Intercept") 
    } else if (threshold == "equidistant") {
      str_add(out$par) <- paste0(
        "  real temp_Intercept1;  // threshold 1 \n",
        "  real", prior[prior$class == "delta", "bound"],
        " delta;  // distance between thresholds \n"
      )
      str_add(out$tparD) <- intercept
      str_add(out$tparC1) <- paste0(
        "  // compute equidistant thresholds \n",
        "  for (k in 1:(ncat - 1)) { \n",
        "    temp_Intercept[k] = temp_Intercept1 + (k - 1.0) * delta; \n",
        "  } \n"
      )
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "temp_Intercept1"), 
        stan_prior(prior, class = "delta")
      )
    }
    
    # generate Stan code specific for each ordinal model
    if (!(family == "cumulative" && ilink == "inv_logit") || disc) {
      cs_arg <- ifelse(!cs, "", "row_vector mucs, ")
      str_add(out$fun) <- paste0(
        "  /* ", family, " log-PDF for a single response \n",
        "   * Args: \n",
        "   *   y: response category \n",
        "   *   mu: linear predictor \n",
        "   *   mucs: optional predictor for category specific effects \n",
        "   *   thres: ordinal thresholds \n",
        "   *   disc: discrimination parameter \n",
        "   * Returns: \n", 
        "   *   a scalar to be added to the log posterior \n",
        "   */ \n",
        "   real ", family, "_lpmf(int y, real mu, ", cs_arg, 
                                  "vector thres, real disc) { \n",
        "     int ncat; \n",
        "     vector[num_elements(thres) + 1] p; \n",
        if (family != "cumulative") 
          "     vector[num_elements(thres)] q; \n",
        "     ncat = num_elements(thres) + 1; \n"
      )
      
      # define actual function content
      if (family == "cumulative") {
        str_add(out$fun) <- paste0(
          "     p[1] = ", ilink, "(", th(1), "); \n",
          "     for (k in 2:(ncat - 1)) { \n", 
          "       p[k] = ", ilink, "(", th("k"), ") - \n",
          "              ", ilink, "(", th("k - 1"), "); \n", 
          "     } \n",
          "     p[ncat] = 1 - ",ilink, "(", th("ncat - 1"), "); \n"
        )
      } else if (family %in% c("sratio", "cratio")) {
        sc <- ifelse(family == "sratio", "1 - ", "")
        str_add(out$fun) <- paste0(
          "     for (k in 1:(ncat - 1)) { \n",
          "       q[k] = ", sc, ilink, "(", th("k"), "); \n",
          "       p[k] = 1 - q[k]; \n",
          "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n", 
          "     } \n",
          "     p[ncat] = prod(q); \n"
        )
      } else if (family == "acat") {
        if (ilink == "inv_logit") {
          str_add(out$fun) <- paste0(
            "     p[1] = 1.0; \n",
            "     for (k in 1:(ncat - 1)) { \n",
            "       q[k] = ", th("k"), "; \n",
            "       p[k + 1] = q[1]; \n",
            "       for (kk in 2:k) p[k + 1] = p[k + 1] + q[kk]; \n",
            "       p[k + 1] = exp(p[k + 1]); \n",
            "     } \n",
            "     p = p / sum(p); \n"
          )
        } else {
          str_add(out$fun) <- paste0(   
            "     for (k in 1:(ncat - 1)) \n",
            "       q[k] = ", ilink, "(", th("k"), "); \n",
            "     for (k in 1:ncat) { \n",     
            "       p[k] = 1.0; \n",
            "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n",
            "       for (kk in k:(ncat - 1)) p[k] = p[k] * (1 - q[kk]); \n",   
            "     } \n",
            "     p = p / sum(p); \n"
          )
        }
      }
      str_add(out$fun) <- "    return categorical_lpmf(y | p); \n   } \n"
    }
  }
  out
}

stan_families <- function(family, bterms) {
  # include .stan files of certain response distributions
  # Args:
  #   family: the model family
  #   bterms: object of class brmsterms
  # Returns:
  #   a list of character strings
  stopifnot(is.family(family), is.brmsterms(bterms))
  families <- family_names(family)
  out <- list()
  if (any(families %in% "categorical")) {
    str_add(out$data) <- "  int<lower=2> ncat;  // number of categories \n" 
    str_add(out$tdataD) <- "  vector[1] zero; \n"
    str_add(out$tdataC) <- "  zero[1] = 0; \n"
  } else if (any(families %in% "zero_inflated_poisson")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_poisson.stan' \n"
  } else if (any(families %in% "zero_inflated_negbinomial")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_negbinomial.stan' \n"
  } else if (any(families %in% "zero_inflated_binomial")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_binomial.stan' \n"
  } else if (any(families %in% "zero_inflated_beta")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_beta.stan' \n"
  } else if (any(families %in% "zero_one_inflated_beta")) {
    str_add(out$fun) <- "  #include 'fun_zero_one_inflated_beta.stan' \n"
  } else if (any(families %in% "hurdle_poisson")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_poisson.stan' \n"
  } else if (any(families %in% "hurdle_negbinomial")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_negbinomial.stan' \n"
  } else if (any(families %in% "hurdle_gamma")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_gamma.stan' \n"
  } else if (any(families %in% "hurdle_lognormal")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_lognormal.stan' \n"
  } else if (any(families %in% "inverse.gaussian")) {
    str_add(out$fun) <- "  #include 'fun_inv_gaussian.stan' \n"
    str_add(out$tdataD) <- "  #include 'tdataD_inv_gaussian.stan' \n"
    str_add(out$tdataC) <- "  #include 'tdataC_inv_gaussian.stan' \n"
  } else if (any(families %in% "von_mises")) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_tan_half.stan' \n",
      "  #include 'fun_von_mises.stan' \n"
    )
  } else if (any(families %in% "wiener")) {
    str_add(out$fun) <- "  #include 'fun_wiener_diffusion.stan' \n"
    str_add(out$tdataD) <- "  real min_Y = min(Y); \n"
  } else if (any(families %in% "asym_laplace")) {
    str_add(out$fun) <- "  #include 'fun_asym_laplace.stan' \n"
  } else if (any(families %in% "skew_normal")) {
    # as suggested by Stephen Martin use sigma and mu of CP 
    # but the skewness parameter alpha of DP
    str_add(out$tdataD) <- "  real sqrt_2_div_pi = sqrt(2 / pi()); \n"
    ap_names <- names(bterms$dpars)
    for (i in which(families %in% "skew_normal")) {
      id <- ifelse(length(families) == 1L, "", i)
      ns <- ifelse(paste0("sigma", id) %in% ap_names, "[n]", "")
      na <- ifelse(paste0("alpha", id) %in% ap_names, "[n]", "")
      type_delta <- ifelse(nzchar(na), "vector[N]", "real")
      no <- ifelse(any(nzchar(c(ns, na))), "[n]", "")
      type_omega <- ifelse(nzchar(no), "vector[N]", "real")
      str_add(out$modelD) <- paste0(
        "  ", type_delta, " delta", id, "; \n",
        "  ", type_omega, " omega", id, "; \n"
      )
      comp_delta <- paste0(
        "  delta", id, na, " = alpha", id, na, 
        " / sqrt(1 + alpha", id, na, "^2); \n"
      )
      comp_omega <- paste0(
        "  omega", id, no, " = sigma", id, ns, 
        " / sqrt(1 - sqrt_2_div_pi^2 * delta", id, na, "^2); \n"
      )
      str_add(out$modelC) <- paste0(
        if (!nzchar(na)) comp_delta,
        if (!nzchar(no)) comp_omega,
        "  for (n in 1:N) { \n",
        if (nzchar(na)) paste0("  ", comp_delta),
        if (nzchar(no)) paste0("  ", comp_omega),
        "    mu", id, "[n] = mu", id, "[n] - omega", id, no, 
        " * delta", id, na, " * sqrt_2_div_pi; \n",
        "  } \n"
      )
    }
  } else if (any(families %in% "gen_extreme_value")) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_gen_extreme_value.stan' \n",
      "  #include 'fun_scale_xi.stan' \n"
    )
    ap_names <- c(names(bterms$dpars), names(bterms$fdpars))
    for (i in which(families %in% "gen_extreme_value")) {
      id <- ifelse(length(families) == 1L, "", i)
      xi <- paste0("xi", id)
      if (!xi %in% ap_names) {
        str_add(out$modelD) <- paste0(
           "  real ", xi, ";  // scaled shape parameter \n"
        )
        sigma <- paste0("sigma", id)
        v <- ifelse(sigma %in% names(bterms$dpars), "_vector", "")
        args <- sargs(paste0("temp_", xi), "Y", paste0("mu", id), sigma)
        str_add(out$modelC) <- paste0(
           "  ", xi, " = scale_xi", v, "(", args, "); \n"
        )
      }
    }
  }
  out
}

stan_mixture <- function(bterms, prior) {
  # Stan code specific for mixture families
  out <- list()
  if (is.mixfamily(bterms$family)) {
    nmix <- length(bterms$family$mix)
    theta_pred <- grepl("^theta", names(bterms$dpars))
    theta_pred <- bterms$dpars[theta_pred]
    theta_fix <- grepl("^theta", names(bterms$fdpars))
    theta_fix <- bterms$fdpars[theta_fix]
    def_thetas <- collapse(
      "  real<lower=0,upper=1> theta", 1:nmix, ";",
      "  // mixing proportion \n"
    )
    if (length(theta_pred)) {
      if (length(theta_pred) != nmix - 1) {
        stop2("Can only predict all but one mixing proportion.")
      }
      missing_id <- setdiff(1:nmix, dpar_id(names(theta_pred)))
      str_add(out$modelD) <- paste0(
        "  vector[N] theta", missing_id, " = rep_vector(0, N); \n",                   
        "  real log_sum_exp_theta; \n"      
      )
      sum_exp_theta <- paste0(
        "exp(theta", 1:nmix, "[n])", collapse = " + "
      )
      str_add(out$modelC3) <- paste0( 
        "    log_sum_exp_theta = log(", sum_exp_theta, "); \n",
        collapse(
          "    theta", 1:nmix, "[n] = theta", 1:nmix, "[n]", 
          " - log_sum_exp_theta; \n"
        )
      )
    } else if (length(theta_fix)) {
      if (length(theta_fix) != nmix) {
        stop2("Can only fix no or all mixing proportions.")
      }
      str_add(out$data) <- paste0(
        "  // mixing proportions \n",                
        collapse(
          "  real<lower=0,upper=1> theta", 1:nmix, "; \n"
        )
      )
    } else {
      str_add(out$data) <- paste0(
        "  vector[", nmix, "] con_theta;  // prior concentration \n"                  
      )
      str_add(out$par) <- paste0(
        "  simplex[", nmix, "] theta;",
        "  // mixing proportions \n"
      )
      str_add(out$prior) <- paste0(
        "  target += dirichlet_lpdf(theta | con_theta); \n"                
      )
      str_add(out$tparD) <- paste0(
        "  // mixing proportions \n",                
        collapse(
          "  real<lower=0,upper=1> theta", 1:nmix, 
          " = theta[", 1:nmix, "]; \n"
        )
      )
    }
    if (bterms$family$order %in% "mu") {
      str_add(out$par) <- paste0( 
        "  ordered[", nmix, "] ordered_Intercept;",
        "  // to identify mixtures \n"
      )
    }
  }
  out
}

stan_se <- function(se) {
  out <- list()
  if (se) {
    str_add(out$data) <- "  vector<lower=0>[N] se;  // known sampling error \n"
    str_add(out$tdataD) <- "  vector<lower=0>[N] se2 = square(se); \n"
  }
  out
}

stan_cens <- function(cens, family) {
  out <- list()
  if (cens) {
    stopifnot(is.family(family))
    str_add(out$data) <- paste0(
      "  int<lower=-1,upper=2> cens[N];  // indicates censoring \n",
      if (isTRUE(attr(cens, "interval"))) {
        paste0(
          ifelse(use_int(family), " int rcens[N];", "  vector[N] rcens;"),
          "  // right censor points for interval censoring \n"
        )
      }
    )
  }
  out
}

stan_disp <- function(bterms, family) {
  # stan code for models with addition argument 'disp'
  # Args:
  #   bterms: object of class brmsterms
  #   family: the model family
  stopifnot(is.brmsterms(bterms))
  stopifnot(is.family(family))
  out <- list()
  if (is.formula(bterms$adforms$disp)) {
    warning2("Addition argument 'disp' is deprecated. ",
             "See help(brmsformula) for more details.")
    par <- if (has_sigma(family)) "sigma"
           else if (has_shape(family)) "shape"
    if (!is.null(bterms[[par]])) {
      stop2("Specifying 'disp' is not allowed when predicting '", par, "'.")
    }
    str_add(out$data) <- "  vector<lower=0>[N] disp;  // dispersion factors \n"
    str_add(out$modelD) <- paste0("  vector[N] disp_", par, "; \n")
    str_add(out$modelC1) <- paste0("  disp_", par, " = ", par, " * disp; \n")
  }
  out
}

stan_pred_functions <- function(x) {
  # add certain predictor functions to Stan's functions block
  out <- ""
  if (grepl("[^[:alnum:]]mo\\(", collapse(x))) {
    str_add(out) <- "  #include fun_monotonic.stan \n"
  } 
  if (grepl("[^[:alnum:]]gp\\(", collapse(x))) {
    str_add(out) <- "  #include fun_gaussian_process.stan \n"
  }
  out
}

stan_misc_functions <- function(family, prior, kronecker) {
  # stan code for user defined functions
  # Args:
  #   family: the model family
  #   prior: object of class brmsprior
  #   kronecker: logical; is the kronecker product needed?
  # Returns:
  #   a string containing defined functions in stan code
  stopifnot(is.family(family))
  out <- ""
  if (family$link == "cauchit") {
    str_add(out) <- "  #include 'fun_cauchit.stan' \n"
  } else if (family$link == "cloglog") {
    str_add(out) <- "  #include 'fun_cloglog.stan' \n"
  }
  if (family$family %in% c("student", "frechet")) {
    str_add(out) <- "  #include 'fun_logm1.stan' \n"
  }
  hs_dfs <- ulapply(attr(prior, "special"), "[[", "hs_df")
  if (any(nzchar(hs_dfs))) {
    str_add(out) <- "  #include 'fun_horseshoe.stan' \n"
  }
  if (kronecker) {
    str_add(out) <- paste0(
      "  #include 'fun_as_matrix.stan' \n",
      "  #include 'fun_kronecker.stan' \n"
    )
  }
  out
}

stan_prior <- function(prior, class, coef = "", group = "", 
                       px = list(), prefix = "", suffix = "", 
                       wsp = 2, matrix = FALSE) {
  # Define priors for parameters in Stan language
  # Args:
  #   prior: an object of class 'brmsprior'
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   nlpar: the name of a non-linear parameter
  #   prefix: a prefix to put at the parameter class
  #   suffix: a suffix to put at the parameter class
  #   matrix: logical; corresponds the class to a parameter matrix?
  #   wsp: an integer >= 0 defining the number of spaces 
  #        in front of the output string
  # Returns:
  #   A character strings in stan language that defines priors 
  #   for a given class of parameters. If a parameter has has 
  #   no corresponding prior in prior, an empty string is returned.
  tp <- tp(wsp)
  wsp <- wsp(nsp = wsp)
  prior_only <- identical(attr(prior, "sample_prior"), "only")
  prior <- subset2(prior, 
    class = class, coef = c(coef, ""), group = c(group, "")
  )
  if (class %in% c("sd", "cor")) {
    # only sd and cor parameters have global priors
    px_tmp <- lapply(px, function(x) c(x, ""))
    prior <- subset2(prior, ls = px_tmp)
  } else {
    prior <- subset2(prior, ls = px)
  }
  if (!nchar(class) && nrow(prior)) {
    # unchecked prior statements are directly passed to Stan
    return(collapse(wsp, prior$prior, "; \n"))
  } 
  
  px <- as.data.frame(px)
  upx <- unique(px)
  if (nrow(upx) > 1L) {
    # can only happen for SD parameters of the same ID
    base_prior <- rep(NA, nrow(upx))
    for (i in seq_len(nrow(upx))) {
      sub_upx <- lapply(upx[i, ], function(x) c(x, ""))
      sub_prior <- subset2(prior, ls = sub_upx) 
      base_prior[i] <- stan_base_prior(sub_prior)
    }
    if (length(unique(base_prior)) > 1L) {
      # define prior for single coefficients manually
      # as there is not single base_prior anymore
      prior_of_coefs <- prior[nzchar(prior$coef), vars_prefix()]
      take <- match_rows(prior_of_coefs, upx)
      prior[nzchar(prior$coef), "prior"] <- base_prior[take]
    }
    base_prior <- base_prior[1]
    bound <- ""
  } else {
    base_prior <- stan_base_prior(prior)
    bound <- prior[!nzchar(prior$coef), "bound"]
  }
  
  individual_prior <- function(i, prior, max_index) {
    # individual priors for each parameter of a class
    if (max_index > 1L || matrix) {
      index <- paste0("[", i, "]")      
    } else {
      index <- ""
    }
    if (nrow(px) > 1L) {
      prior <- subset2(prior, ls = px[i, ])
    }
    uc_prior <- prior$prior[match(coef[i], prior$coef)]
    if (!is.na(uc_prior) & nchar(uc_prior)) { 
      # user defined prior for this parameter
      coef_prior <- uc_prior
    } else { 
      # base prior for this parameter
      coef_prior <- base_prior 
    }  
    if (nzchar(coef_prior)) {  
      # implies a proper prior
      pars <- paste0(class, index)
      out <- stan_target_prior(coef_prior, pars, bound = bound)
      out <- paste0(tp, out, "; \n")
    } else {
      # implies an improper flat prior
      out <- ""
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(prefix, class, suffix)
  if (any(with(prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(
      seq_along(coef), individual_prior, 
      prior = prior, max_index = length(coef)
    )
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- paste0("to_vector(", class, ")")
    }
    out <- stan_target_prior(
      base_prior, class, ncoef = length(coef), bound = bound
    )
    out <- paste0(tp, out, "; \n")
  } else {
    out <- ""
  }
  special_prior <- stan_special_prior(
    class, prior, ncoef = length(coef), px = px
  )
  out <- collapse(c(out, special_prior))
  if (prior_only && nzchar(class) && !nchar(out)) {
    stop2("Sampling from priors is not possible as ", 
          "some parameters have no proper priors. ",
          "Error occured for class '", class, "'.")
  }
  out
}

stan_base_prior <- function(prior) {
  # get base (highest level) prior of all priors
  # Args:
  #   prior: a prior.frame
  stopifnot(length(unique(prior$class)) <= 1L) 
  prior <- prior[with(prior, !nzchar(coef) & nzchar(prior)), ]
  vars <- c("group", "nlpar", "dpar", "resp", "class")
  i <- 1
  found <- FALSE
  base_prior <- ""
  take <- rep(FALSE, nrow(prior))
  while (!found && i <= length(vars)) {
    take <- nzchar(prior[[vars[i]]]) & !take
    if (any(take)) {
      base_prior <- prior[take, "prior"]
      found <- TRUE
    }
    i <- i + 1
  }
  stopifnot(length(base_prior) == 1L)
  base_prior
}

stan_target_prior <- function(prior, par, ncoef = 1, bound = "") {
  prior_name <- get_matches("^[^\\(]+\\(", prior, simplify = FALSE)
  for (i in seq_along(prior_name)) {
    if (length(prior_name[[i]]) != 1L) {
      stop2("The prior '", prior[i], "' is invalid.")
    }
  }
  prior_name <- unlist(prior_name)
  prior_name <- substr(prior_name, 1, nchar(prior_name) - 1)
  prior_args <- rep(NA, length(prior))
  for (i in seq_along(prior)) {
    prior_args[i] <- sub(
      paste0("^", prior_name[i], "\\("), "", prior[i]
    )
  }
  out <- paste0(prior_name, "_lpdf(", par, " | ", prior_args)
  par_class <- unique(get_matches("^[^_]+", par))
  par_bound <- par_bounds(par_class, bound)
  prior_bound <- prior_bounds(prior_name)
  trunc_lb <- is.character(par_bound$lb) || par_bound$lb > prior_bound$lb
  trunc_ub <- is.character(par_bound$ub) || par_bound$ub < prior_bound$ub
  if (trunc_lb || trunc_ub) {
    wsp <- wsp(nsp = 4)
    if (trunc_lb && !trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lccdf(", 
        par_bound$lb, " | ", prior_args
      )
    } else if (!trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lcdf(", 
        par_bound$ub, " | ", prior_args
      )
    } else if (trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * log_diff_exp(", 
        prior_name, "_lcdf(", par_bound$ub, " | ", prior_args, ", ",
        prior_name, "_lcdf(", par_bound$lb, " | ", prior_args, ")"
      )
    }
  }
  out
}

stan_special_prior <- function(class, prior, ncoef, px = list()) {
  # add special priors such as horseshoe and lasso
  out <- ""
  p <- usc(combine_prefix(px))
  if (all(class == paste0("b", p))) {
    stopifnot(length(p) == 1L)
    tp <- tp()
    # add horseshoe and lasso shrinkage priors
    prefix <- combine_prefix(px, keep_mu = TRUE)
    special <- attributes(prior)$special[[prefix]]
    if (!is.null(special$hs_df)) {
      local_args <- paste0("0.5 * hs_df", p)
      local_args <- sargs(local_args, local_args)
      global_args <- paste0("0.5 * hs_df_global", p)
      global_args <- sargs(global_args, global_args)
      c2_args <- paste0("0.5 * hs_df_slab", p)
      c2_args <- sargs(c2_args, c2_args)
      wsp <- wsp(nsp = 4)
      str_add(out) <- paste0(
        tp, "normal_lpdf(zb", p, " | 0, 1); \n",
        tp, "normal_lpdf(hs_local", p, "[1] | 0, 1)\n", 
        wsp, "- ", ncoef, " * log(0.5); \n",
        tp, "inv_gamma_lpdf(hs_local", p, "[2] | ", local_args, "); \n",
        tp, "normal_lpdf(hs_global", p, "[1] | 0, 1)\n", 
        wsp, "- 1 * log(0.5); \n",
        tp, "inv_gamma_lpdf(hs_global", p, "[2] | ", global_args, "); \n",
        tp, "inv_gamma_lpdf(hs_c2", p, " | ", c2_args, "); \n"
      )
    }
    if (!is.null(special$lasso_df)) {
      str_add(out) <- paste0(
        tp, "chi_square_lpdf(lasso_inv_lambda", p, " | lasso_df", p, "); \n"
      )
    }
  }
  out
}

stan_rngprior <- function(sample_prior, prior, par_declars,
                          gen_quantities, prior_special) {
  # stan code to sample from priors seperately
  # Args:
  #   sample_prior: take samples from priors?
  #   prior: character string taken from stan_prior
  #   par_declars: the parameters block of the Stan code
  #     required to extract boundaries
  #   gen_quantities: Stan code from the generated quantities block
  #   prior_special: a list of values pertaining to special priors
  #     such as horseshoe or lasso
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  out <- list()
  if (identical(sample_prior, "yes")) {
    prior <- strsplit(gsub(" |\\n", "", prior), ";")[[1]]
    prior <- prior[nzchar(prior)]
    pars <- get_matches("_lpdf\\([^|]+", prior, first = TRUE)
    pars <- gsub("^_lpdf\\(|to_vector\\(|\\)$", "", pars)
    regex <- "^(z|zs|zb|zgp|Xme|hs)_?|^increment_log_prob\\("
    take <- !grepl(regex, pars)
    pars <- rename(
      pars[take], symbols = c("^L_", "^Lrescor"), 
      subs = c("cor_", "rescor"), fixed = FALSE
    )
    dis <- get_matches("target\\+=[^\\(]+", prior, first = TRUE)
    dis <- gsub("target\\+=|_lpdf", "", dis[take])
    args <- get_matches("\\|[^$\\|]+\\)($|-)", prior, first = TRUE)
    args <- gsub("\\||(;|-)$", "", args[take])
    type <- rep("real", length(pars))
    
    # rename parameters containing indices
    has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
    pars[has_ind] <- ulapply(pars[has_ind], function(par) {
      ind <- regmatches(par, gregexpr("\\[[[:digit:]]+\\]", par))
      ind <- as.numeric(substr(ind, 2, nchar(ind) - 1))
      gsub("\\[[[:digit:]]+\\]", paste0("_", ind), par)
    })
    
    # special treatment of lkj_corr_cholesky priors
    args <- ifelse(
      grepl("corr_cholesky$", dis), 
      paste0("2,", args, "[1, 2]"), args
    )
    dis <- sub("corr_cholesky$", "corr", dis)
    
    # extract information from the initial parameter definition
    par_declars <- unlist(strsplit(par_declars, "\n", fixed = TRUE))
    par_declars <- gsub("^[[:blank:]]*", "", par_declars)
    par_declars <- par_declars[!grepl("^//", par_declars)]
    all_pars <- get_matches(" [^[:blank:]]+;", par_declars) 
    all_pars <- substr(all_pars, 2, nchar(all_pars) - 1)
    all_bounds <- get_matches("<.+>", par_declars, simplify = FALSE)
    all_bounds <- ulapply(all_bounds, function(x) if (length(x)) x else "")
    all_types <- get_matches("^[^[:blank:]]+", par_declars)
    
    # define parameter types and boundaries
    bounds <- rep("", length(pars))
    types <- rep("real", length(pars))
    for (i in seq_along(all_pars)) {
      k <- which(grepl(paste0("^", all_pars[i]), pars))
      bounds[k] <- all_bounds[i]
      if (grepl("^simplex", all_pars[i])) {
        types[k] <- all_types[i]
      }
    }
    
    # distinguish between bounded and unbounded parameters
    has_bounds <- as.logical(nchar(bounds))
    if (any(has_bounds)) {  
      # bounded parameters have to be sampled in the model block
      str_add(out$par) <- paste0(
        "  // parameters to store prior samples\n",
        collapse(
          "  real", bounds[has_bounds], 
          " prior_", pars[has_bounds], ";\n"
        )
      )
      str_add(out$model) <- paste0(
        "  // additionally draw samples from priors\n",
        collapse(
          "  target += ", dis[has_bounds], "_lpdf(",
          "prior_", pars[has_bounds], " | ", args[has_bounds], "; \n"
        )
      )
    }
    no_bounds <- !has_bounds
    if (any(no_bounds)) {
      # use parameters sampled from priors for use in other priors
      spars <- NULL
      # cannot sample from the horseshoe prior anymore as of brms 1.5.0
      lasso_prefix <- nzchar(ulapply(prior_special, "[[", "lasso_df"))
      lasso_prefix <- names(prior_special)[lasso_prefix]
      lasso_prefix <- usc(sub("^mu(_|$)", "", lasso_prefix))
      if (length(lasso_prefix)) {
        spars <- c(spars, paste0("lasso_inv_lambda", lasso_prefix))
      }
      if (length(spars)) {
        bpars <- grepl("^b(|mo|cs|me)(_|$)", pars)
        args[bpars] <- rename(args[bpars], spars, paste0("prior_", spars))
      }
      # unbounded parameters can be sampled in the generatated quantities block
      str_add(out$genD) <- paste0(
        "  // additionally draw samples from priors \n",
        collapse(
          "  ", types[no_bounds], " prior_", pars[no_bounds], 
          " = ", dis[no_bounds], "_rng(", args[no_bounds], "; \n"
        )
      )
    }
    # compute priors for the actual population-level intercepts
    is_temp_intercept <- grepl("^temp.*_Intercept", pars)
    if (any(is_temp_intercept)) {
      temp_intercepts <- pars[is_temp_intercept]
      p <- gsub("^temp|_Intercept$", "", temp_intercepts)
      intercepts <- paste0("b", p, "_Intercept")
      regex <- paste0(" (\\+|-) dot_product\\(means_X", p, ", b", p, "\\)")
      sub_X_means <- lapply(regex, get_matches, gen_quantities)
      sub_X_means <- ulapply(sub_X_means, 
        function(x) if (length(x)) x[1] else ""
      )
      str_add(out$genD) <- collapse(
        "  real prior_", intercepts,
        " = prior_", temp_intercepts, sub_X_means, "; \n"
      )
    }
  }
  out
}

stan_link <- function(link) {
  # find the link in Stan language
  # Args:
  #   link: the link function
  switch(link, 
    identity = "",
    log = "log", 
    logm1 = "logm1",
    inverse = "inv",
    sqrt = "sqrt", 
    "1/mu^2" = "inv_square", 
    logit = "logit", 
    probit = "inv_Phi", 
    probit_approx = "inv_Phi", 
    cloglog = "cloglog", 
    cauchit = "cauchit",
    tan_half = "tan_half",
    log1p = "log1p"
  )
}

stan_ilink <- function(link) {
  # find the inverse link in Stan language
  # Args:
  #   link: the link function
  switch(link, 
    identity = "",
    log = "exp", 
    logm1 = "expp1",
    inverse = "inv", 
    sqrt = "square", 
    "1/mu^2" = "inv_sqrt", 
    logit = "inv_logit", 
    probit = "Phi", 
    probit_approx = "Phi_approx", 
    cloglog = "inv_cloglog",
    cauchit = "inv_cauchit",
    tan_half = "inv_tan_half",
    log1p = "expm1"
  )
}

stan_has_built_in_fun <- function(family) {
  # indicates if a family-link combination has a build in 
  # function in Stan (such as binomial_logit)
  # Args:
  #   family: a list with elements 'family' and 'link'
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family$link
  par <- family$par
  family <- family$family
  log_families <- c(
    "poisson", "negbinomial", "geometric", 
    "zero_inflated_poisson", "zero_inflated_negbinomial",
    "hurdle_poisson", "hurdle_negbinomial"
  )
  logit_families <- c(
    "binomial", "bernoulli", "cumulative", "categorical",
    "zero_inflated_binomial"
  )
  isTRUE(
    family %in% log_families && link == "log" ||
    family %in% logit_families && link == "logit" ||
    isTRUE(par %in% c("zi", "hu")) && link == "logit"
  )
}

stan_needs_kronecker <- function(ranef, names_cov_ranef) {
  # checks if a model needs the kronecker product
  # Args: 
  #   ranef: named list returned by tidy_ranef
  #   names_cov_ranef: names of the grouping factors that
  #                    have a cov.ranef matrix 
  ids <- unique(ranef$id)
  out <- FALSE
  for (id in ids) {
    r <- ranef[ranef$id == id, ]
    out <- out || nrow(r) > 1L && r$cor[1] && 
             r$group[1] %in% names_cov_ranef
  }
  out
}
