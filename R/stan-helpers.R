stan_llh <- function(family, effects = list(), data = NULL, 
                     autocor = cor_arma()) {
  # Likelihood in Stan language
  # Args:
  #   family: the model family
  #   effects: output of extract_effects
  #   autocor: object of classe cor_brms
  stopifnot(is(family, "family"))
  link <- family$link
  family <- family$family
  is_categorical <- is.categorical(family)
  is_ordinal <- is.ordinal(family)
  is_hurdle <- is.hurdle(family)
  is_zero_inflated <- is.zero_inflated(family)
  is_forked <- is.forked(family)
  is_mv <- is.linear(family) && length(effects$response) > 1L
  
  has_se <- is.formula(effects$se)
  has_weights <- is.formula(effects$weights)
  has_cens <- has_cens(effects$cens, data = data)
  has_disp <- is.formula(effects$disp)
  has_trials <- is.formula(effects$trials)
  has_cse <- has_cse(effects)
  bounds <- get_bounds(effects$trunc, data = data)
  has_trunc <- any(bounds$lb > -Inf) || any(bounds$ub < Inf)
  ll_adj <- has_cens || has_weights || has_trunc

  if (is_mv) {
    # prepare for use of a multivariate likelihood
    family <- paste0(family, "_mv")
  } else if (use_cov(autocor)) {
    # ARMA effects in covariance matrix formulation
    if (ll_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    family <- paste0(family, "_cov")
  } else if (is(autocor, "cov_fixed")) {
    if (has_se || ll_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    family <- paste0(family, "_fixed")
  }
  
  auxpars <- intersect(auxpars(), names(effects))
  reqn <- ll_adj || is_categorical || is_ordinal || 
          is_hurdle || is_zero_inflated || 
          is.exgaussian(family) || is.wiener(family) ||
          any(c("phi", "kappa") %in% auxpars)
  n <- ifelse(reqn, "[n]", "")
  # prepare auxiliary parameters
  disp <- ifelse(has_disp, "disp_", "")
  reqn_sigma <- (ll_adj && (has_se || has_disp)) || 
                (reqn && "sigma" %in% auxpars)
  sigma <- paste0(ifelse(has_se, "se", paste0(disp, "sigma")),
                  if (reqn_sigma) "[n]")
  reqn_shape <- (ll_adj && has_disp) || (reqn && "shape" %in% auxpars)
  shape <- paste0(disp, "shape", if (reqn_shape) "[n]")
  for (ap in setdiff(auxpars(), c("sigma", "shape"))) {
    assign(ap, stan_apn(ap, auxpars, reqn))
  }
  .logit <- ifelse(any(c("zi", "hu") %in% auxpars), "_logit", "")
  reqn_trials <- has_trials && (ll_adj || is_zero_inflated)
  trials <- ifelse(reqn_trials, "trials[n]", "trials")

  simplify <- !has_trunc && !has_cens && stan_has_built_in_fun(family, link)
  eta <- paste0(ifelse(is_mv, "Eta", "eta"), n)
  ord_args <- paste0("eta[n], ", if (has_cse) "etacs[n], ", "temp_Intercept")
  inv_gauss_fun <- paste0("inv_gaussian", if (!reqn) "_vector")
  inv_gauss_args <- paste0(eta, ", shape, ", if (!reqn) "sum_", 
                           "log_Y", n, ", sqrt_Y", n)
  if (simplify) { 
    llh_pre <- switch(family,
      poisson = c("poisson_log", eta), 
      negbinomial = c("neg_binomial_2_log", paste0(eta, ", ", shape)),
      geometric = c("neg_binomial_2_log", paste0(eta, ", 1")),
      cumulative = c("ordered_logistic", "eta[n], temp_Intercept"),
      categorical = c("categorical_logit", "append_row(zero, Eta[n])"),
      binomial = c("binomial_logit", paste0(trials, ", ", eta)), 
      bernoulli = c("bernoulli_logit", eta))
  } else {
    llh_pre <- switch(family,
      gaussian = c("normal", paste0(eta, ", ", sigma)),
      gaussian_cov = c("normal_cov", paste0(eta,", se2, N_tg, ", 
                       "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      gaussian_mv = c("multi_normal_cholesky", paste0(eta, ", LSigma")),
      gaussian_fixed = c("multi_normal_cholesky", paste0(eta, ", LV")),
      student = c("student_t",  paste0(nu, ", ", eta, ", ", sigma)),
      student_cov = c("student_t_cov", paste0(nu, ", ", eta, ", se2, N_tg, ", 
                      "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      student_mv = c("multi_student_t", paste0(nu, ", ", eta, ", Sigma")),
      student_fixed = c("multi_student_t", paste0(nu, ", ", eta, ", V")),
      lognormal = c("lognormal", paste0(eta, ", ", sigma)),
      poisson = c("poisson", eta),
      negbinomial = c("neg_binomial_2", paste0(eta, ", ", shape)),
      geometric = c("neg_binomial_2", paste0(eta, ", 1")),
      binomial = c("binomial", paste0(trials, ", ", eta)),
      bernoulli = c("bernoulli", eta),
      gamma = c("gamma", paste0(shape, ", ", eta)), 
      exponential = c("exponential", eta),
      weibull = c("weibull", paste0(shape, ", ", eta)), 
      exgaussian = c("exgaussian", paste0(eta, ", ", sigma, ", ", beta)),
      inverse.gaussian = c(inv_gauss_fun, inv_gauss_args),
      wiener <- c("wiener_diffusion", 
                  paste0("dec[n], ", eta, ", ", bs, ", ", ndt, ", ", bias)),
      beta = c("beta", paste0(eta, " * ", phi, ", (1 - ", eta, ") * ", phi)),
      von_mises = c(paste0("von_mises_", ifelse(reqn, "real", "vector")), 
                           paste0(eta, ", ", kappa)),
      cumulative = c("cumulative", ord_args),
      sratio = c("sratio", ord_args),
      cratio = c("cratio", ord_args),
      acat = c("acat", ord_args),
      hurdle_poisson = c(paste0("hurdle_poisson", .logit), 
                         paste0(eta, ", ", hu)),
      hurdle_negbinomial = c(paste0("hurdle_neg_binomial", .logit), 
                             paste0(eta, ", ", hu, ", ", shape)),
      hurdle_gamma = c(paste0("hurdle_gamma", .logit), 
                       paste0(shape, ", ", eta, ", ", hu)),
      hurdle_lognormal = c(paste0("hurdle_lognormal", .logit), 
                           paste0(eta, ", ", hu, ", ", sigma)),
      zero_inflated_poisson = c(paste0("zero_inflated_poisson", .logit), 
                                paste0(eta, ", ", zi)),
      zero_inflated_negbinomial = 
        c(paste0("zero_inflated_neg_binomial", .logit),
          paste0(eta, ", ", zi, ", ", shape)),
      zero_inflated_binomial = c(paste0("zero_inflated_binomial", .logit), 
                                 paste0(trials, ", ", eta, ", ", zi)),
      zero_inflated_beta = c(paste0("zero_inflated_beta", .logit), 
                             paste0(eta, ", ", zi, ", ", phi)))
  }
  
  # write likelihood code
  interval <- isTRUE(attr(has_cens, "interval"))
  type <- c("cens", "weights")[match(TRUE, c(has_cens, has_weights))]
  if (is.na(type)) type <- "general"
  llh <- switch(type, 
    cens = stan_llh_cens(llh_pre, family, interval, has_weights),
    weights = stan_llh_weights(llh_pre, family),
    general = stan_llh_general(llh_pre, reqn, bounds)) 
  if (reqn) {
    # loop over likelihood if it cannot be vectorized
    llh <- paste0("  for (n in 1:N) { \n    ", llh, "    } \n")
  }
  llh
}

stan_llh_general <- function(llh_pre, reqn = FALSE, bounds = NULL) {
  # default likelihood in Stan language
  # Args:
  #   reqn: does Y require the index 'n'?
  #   bounds: a list containing elements lb and ub
  stopifnot(length(llh_pre) == 2L)
  if (any(bounds$lb > -Inf) || any(bounds$ub < Inf)) {
    # prepare possible truncation
    lb <- ifelse(any(bounds$lb > -Inf), "lb[n]", "")
    ub <- ifelse(any(bounds$ub < Inf), "ub[n]", "")
    trunc <- paste0(" T[", lb, ", ", ub, "]")
  } else {
    trunc <- ""
  }
  paste0("  Y", ifelse(reqn, "[n]", ""), " ~ ", llh_pre[1], 
         "(", llh_pre[2], ")", trunc, "; \n")
}

stan_llh_cens <- function(llh_pre, family = gaussian(), 
                          interval = FALSE, weights = FALSE) {
  # censored likelihood in Stan language
  # Args: 
  #   interval: are there interval censored responses present?
  #   weights: is the model additionally weighted?
  stopifnot(length(llh_pre) == 2L)
  s <- collapse(rep(" ", 6))
  tp <- "  target += "
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  w <- ifelse(weights, "weights[n] * ", "")
  paste0("  // special treatment of censored data \n",
    s, "if (cens[n] == 0) \n", 
    s, tp, w, llh_pre[1], "_", lpdf, "(Y[n] | ", llh_pre[2], "); \n",
    s, "else if (cens[n] == 1) \n",         
    s, tp, w, llh_pre[1], "_lccdf(Y[n] | ", llh_pre[2], "); \n",
    s, "else if (cens[n] == -1) \n",
    s, tp, w, llh_pre[1], "_lcdf(Y[n] | ", llh_pre[2], "); \n",
    if (interval) {
      paste0(
        s, "else if (cens[n] == 2) \n",
        s, tp, w, "log_diff_exp(", 
          llh_pre[1], "_lcdf(rcens[n] | ", llh_pre[2], "), \n",
        collapse(rep(" ", 31)),
          llh_pre[1], "_lcdf(Y[n] | ", llh_pre[2], ")); \n")
    })
}

stan_llh_weights <- function(llh_pre, family = gaussian()) {
  # weighted likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  paste0("  lp_pre[n] = ", llh_pre[1], "_", lpdf, 
         "(Y[n] | ", llh_pre[2],"); \n")
}

stan_autocor <- function(autocor, effects = list(), family = gaussian(),
                         prior = brmsprior()) {
  # Stan code related to autocorrelation structures
  # Args:
  #   autocor: autocorrelation structure; object of class cor_brms
  #   effects: output of extract_effects
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  stopifnot(is(family, "family"))
  is_linear <- is.linear(family)
  resp <- effects$response
  is_mv <- is_linear && length(resp) > 1L
  link <- stan_link(family$link)
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  Karr <- get_arr(autocor)
  out <- list()
  if (Kar || Kma) {
    if (!is_linear) {
      stop2("The ARMA correlation structure is not yet implemented ", 
            "for family '", family$family, "'.") 
    }
    out$data <- paste0(out$data, "  #include 'data_arma.stan' \n")
    # restrict ARMA effects to be in [-1,1] when using covariance
    # formulation as they cannot be outside this interval anyway
    if (Kar) {
      ar_bound <- with(prior, bound[class == "ar"])
      out$par <- paste0(out$par, 
        "  vector", ar_bound, "[Kar] ar;  // autoregressive effects \n")
      out$prior <- paste0(out$prior, stan_prior(class = "ar", prior = prior))
    }
    if (Kma) {
      ma_bound <- with(prior, bound[class == "ma"])
      out$par <- paste0(out$par, 
        "  vector", ma_bound, "[Kma] ma;  // moving-average effects \n")
      out$prior <- paste0(out$prior, stan_prior(class = "ma", prior = prior))
    }
    
    if (use_cov(autocor)) {
      # if the user wants ARMA effects to be estimated using
      # a covariance matrix for residuals
      err_msg <- "ARMA covariance matrices are not yet working"
      if (is_mv) {
        stop2(err_msg, " in multivariate models.")
      }
      if (is.formula(effects$disp)) {
        stop2(err_msg, " when specifying 'disp'.")
      }
      if ("sigma" %in% names(effects)) {
        stop2(err_msg, " when predicting 'sigma'.")
      }
      out$data <- paste0(out$data, "  #include 'data_arma_cov.stan' \n")
      out$transD <- "  matrix[max(nobs_tg), max(nobs_tg)] res_cov_matrix; \n"
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
      out$transC1 <- paste0("  // compute residual covariance matrix \n",
                            "  res_cov_matrix = cov_matrix_", cov_mat_fun, 
                            "(", cov_mat_args, ", sigma, max(nobs_tg)); \n")
      # defined selfmade functions for the functions block
      if (family$family == "gaussian") {
        out$fun <- paste0(out$fun, "  #include 'fun_normal_cov.stan' \n")
      } else {  # family == "student"
        out$fun <- paste0(out$fun, "  #include 'fun_student_t_cov.stan' \n")
      }
      if (Kar && !Kma) {
        out$fun <- paste0(out$fun, "  #include 'fun_cov_matrix_ar1.stan' \n")
      } else if (!Kar && Kma) {
        out$fun <- paste0(out$fun, "  #include 'fun_cov_matrix_ma1.stan' \n")
      } else {
        out$fun <- paste0(out$fun, "  #include 'fun_cov_matrix_arma1.stan' \n")
      }
    } else {
      err_msg <- "Please set cov = TRUE in cor_arma / cor_ar / cor_ma"
      if (is.formula(effects$se)) {
        stop2(err_msg, " when specifying 'se'.")
      }
      if (length(effects$nonlinear)) {
        stop2(err_msg, " for non-linear models.")
      }
      if (is_mv) {
        rs <- usc(resp, "prefix")
        index <- paste0("n, ", seq_along(rs))
      } else {
        rs <- ""
        index <- "n"
      }
      out$modelD <- paste0(
        "  // objects storing residuals \n",
        collapse("  matrix[N, Karma] E", rs, "; \n",
                 "  vector[N] e", rs, "; \n"))
      out$modelC1 <- collapse("  E", rs, " = rep_matrix(0.0, N, Karma); \n")
      out$modelC2 <- paste0(
        "    // computation of ARMA effects \n",
        collapse("    e", rs, "[n] = ", link, "(Y[", index, "])", 
                 " - eta", rs, "[n]", "; \n"),
        "    for (i in 1:Karma) { \n", 
        "      if (n + 1 - i > 0 && n < N && tg[n + 1] == tg[n + 1 - i]) { \n",
        collapse("        E", rs, "[n + 1, i] = e", rs, "[n + 1 - i]; \n"),
        "      } \n",
        "    } \n")
    } 
  }
  if (Karr) {
    # autoregressive effects of the response
    out$data <- paste0(out$data,
      "  // data needed for ARR effects \n",
      "  int<lower=1> Karr; \n",
      "  matrix[N, Karr] Yarr;  // ARR design matrix \n")
    out$par <- paste0(out$par,
      "  vector", with(prior, bound[class == "arr"]), "[Karr] arr;",
      "  // autoregressive effects of the response \n")
    out$prior <- paste0(out$prior, stan_prior(class = "arr", prior = prior))
  }
  if (is(autocor, "cov_fixed")) {
    if (!is_linear) {
      stop2("Fixed residual covariance matrices are not yet ", 
            "implemented for family '", family$family, "'.") 
    }
    out$data <- "  matrix[N, N] V;  // known residual covariance matrix \n"
    if (family$family %in% "gaussian") {
      out$tdataD <- "  matrix[N, N] LV; \n"
      out$tdataC <- "  LV = cholesky_decompose(V); \n"
    }
  }
  if (is(autocor, "cor_bsts")) {
    if (is_mv || family$family %in% c("bernoulli", "categorical")) {
      stop2("The bsts structure is not yet implemented for this family.")
    }
    if (length(effects$nonlinear)) {
      stop2("The bsts structure is not yet implemented for non-linear models.")
    }
    out$data <- "  vector[N] tg;  // indicates independent groups \n"
    out$par <- paste0("  vector[N] loclev;  // local level terms \n",
                      "  real<lower=0> sigmaLL;  // SD of local level terms \n")
    if (is_linear && !is_mv) {
      # this often helps with convergence
      center <- paste0(link, "(Y[", c("1", "n"), "])")
    } else {
      center <- c("0", "0")
    }
    out$prior <- paste0(out$prior, 
      stan_prior(class = "sigmaLL", prior = prior),
      "  loclev[1] ~ normal(", center[1], ", sigmaLL); \n",
      "  for (n in 2:N) { \n",
      "    if (tg[n] == tg[n - 1]) { \n",
      "      loclev[n] ~ normal(loclev[n - 1], sigmaLL); \n",
      "    } else { \n",
      "      loclev[n] ~ normal(", center[2], ", sigmaLL); \n",
      "    } \n",
      "  } \n")
  }
  out
}

stan_mv <- function(family, response, prior = brmsprior()) {
  # some Stan code for multivariate models
  # Args:
  #   family: model family
  #   response: names of the response variables
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # Returns: 
  #   list containing Stan code specific for multivariate models
  stopifnot(is(family, "family"))
  out <- list()
  nresp <- length(response)
  if (nresp > 1L) {
    if (is.linear(family)) {
      out$data <- "  #include 'data_mv.stan' \n"
      out$par <- paste0(
        "  // parameters for multivariate linear models \n",
        "  vector<lower=0>[nresp] sigma; \n",
        "  cholesky_factor_corr[nresp] Lrescor; \n")
      out$prior <- paste0(
        stan_prior(class = "sigma", coef = response, prior = prior),
        stan_prior(class = "Lrescor", prior = prior))
      if (family$family == "gaussian") {
        out$transD <- "  cholesky_factor_cov[nresp] LSigma; \n"
        out$transC1 <- paste0(
          "  // compute cholesky factor of residual covariance matrix \n",
          "  LSigma = diag_pre_multiply(sigma, Lrescor); \n")
      } else if (family$family == "student") {
        out$transD <- "  cov_matrix[nresp] Sigma; \n"
        out$transC1 <- paste0(
          "  // compute residual covariance matrix \n",
          "  Sigma = multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor)); \n")
      }
      out$genD <- paste0(
        "  matrix[nresp, nresp] Rescor; \n",
        "  vector<lower=-1,upper=1>[nrescor] rescor; \n")
      out$genC <- paste0(
        "  // take only relevant parts of residual correlation matrix \n",
        "  Rescor = multiply_lower_tri_self_transpose(Lrescor); \n",
        collapse(ulapply(2:nresp, function(i) lapply(1:(i-1), function(j)
          paste0("  rescor[",(i-1)*(i-2)/2+j,"] = Rescor[",j,", ",i,"]; \n")))))
    } else if (!is.categorical(family)) {
      stop2("Multivariate models are not yet implemented ", 
            "for family '", family$family, "'.")
    }
  }
  out
}

stan_ordinal <- function(family, prior = brmsprior(), 
                         cse = FALSE, threshold = "flexible") {
  # Ordinal effects in Stan
  # Args:
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   cse: logical; are there category specific effects?
  #   threshold: either "flexible" or "equidistant" 
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  stopifnot(is(family, "family"))
  out <- list()
  if (is.ordinal(family)) {
    # define Stan code similar for all ordinal models
    out$data <- "  int ncat;  // number of categories \n"
    th <- function(k, fam = family) {
      # helper function generating stan code inside ilink(.)
      sign <- ifelse(fam %in% c("cumulative", "sratio"), " - ", " + ")
      ptl <- ifelse(cse, paste0(sign, "etacs[k]"), "") 
      if (sign == " - ") {
        out <- paste0("thres[",k,"]", ptl, " - eta")
      } else {
        out <- paste0("eta", ptl, " - thres[",k,"]")
      }
    } 
    link <- family$link
    family <- family$family
    ilink <- stan_ilink(link)
    type <- ifelse(family == "cumulative", "ordered", "vector")
    intercept <- paste0("  ", type, "[ncat-1] temp_Intercept;",
                        "  // temporary thresholds \n")
    if (threshold == "flexible") {
      out$par <- intercept
      out$prior <- stan_prior("temp_Intercept", prior = prior) 
    } else if (threshold == "equidistant") {
      out$par <- paste0("  real temp_Intercept1;  // threshold 1 \n",
                        "  real", if (family == "cumulative") "<lower=0>",
                        " delta;  // distance between thresholds \n")
      out$transD <- intercept
      out$transC1 <- paste0(
        "  // compute equidistant thresholds \n",
        "  for (k in 1:(ncat - 1)) { \n",
        "    temp_Intercept[k] = temp_Intercept1 + (k - 1.0) * delta; \n",
        "  } \n")
      out$prior <- paste0(stan_prior(class = "temp_Intercept1", prior = prior), 
                          stan_prior(class = "delta", prior = prior))
    }
    
    # generate Stan code specific for each ordinal model
    if (!(family == "cumulative" && ilink == "inv_logit")) {
      cse_arg <- ifelse(!cse, "", "row_vector etacs, ")
      out$fun <- paste0(
        "  /* ", family, " log-PDF for a single response \n",
        "   * Args: \n",
        "   *   y: response category \n",
        "   *   eta: linear predictor \n",
        "   *   etacs: optional predictor for category specific effects \n",
        "   *   thres: ordinal thresholds \n",
        "   * Returns: \n", 
        "   *   a scalar to be added to the log posterior \n",
        "   */ \n",
        "   real ", family, "_lpmf(int y, real eta, ", cse_arg, "vector thres) { \n",
        "     int ncat; \n",
        "     vector[num_elements(thres) + 1] p; \n",
        if (family != "cumulative") "     vector[num_elements(thres)] q; \n",
        "     ncat = num_elements(thres) + 1; \n")
      
      # define actual function content
      if (family == "cumulative") {
        out$fun <- paste0(out$fun,
        "     p[1] = ", ilink, "(", th(1), "); \n",
        "     for (k in 2:(ncat - 1)) { \n", 
        "       p[k] = ", ilink, "(", th("k"), ") - ",
        ilink, "(", th("k - 1"), "); \n", 
        "     } \n",
        "     p[ncat] = 1 - ",ilink, "(", th("ncat - 1"), "); \n")
      } else if (family %in% c("sratio", "cratio")) {
        sc <- ifelse(family == "sratio", "1 - ", "")
        out$fun <- paste0(out$fun,
        "     for (k in 1:(ncat - 1)) { \n",
        "       q[k] = ", sc, ilink, "(", th("k"), "); \n",
        "       p[k] = 1 - q[k]; \n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n", 
        "     } \n",
        "     p[ncat] = prod(q); \n")
      } else if (family == "acat") {
        if (ilink == "inv_logit") {
          out$fun <- paste0(out$fun,
          "     p[1] = 1.0; \n",
          "     for (k in 1:(ncat - 1)) { \n",
          "       q[k] = ", th("k"), "; \n",
          "       p[k + 1] = q[1]; \n",
          "       for (kk in 2:k) p[k + 1] = p[k + 1] + q[kk]; \n",
          "       p[k + 1] = exp(p[k + 1]); \n",
          "     } \n",
          "     p = p / sum(p); \n")
        } else {
          out$fun <- paste0(out$fun,    
          "     for (k in 1:(ncat - 1)) \n",
          "       q[k] = ", ilink, "(", th("k"), "); \n",
          "     for (k in 1:ncat) { \n",     
          "       p[k] = 1.0; \n",
          "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n",
          "       for (kk in k:(ncat - 1)) p[k] = p[k] * (1 - q[kk]); \n",      
          "     } \n",
          "     p = p / sum(p); \n")
        }
      }
      out$fun <- paste(out$fun, 
        "    return categorical_lpmf(y | p); \n   } \n")
    }
  }
  out
}

stan_families <- function(family) {
  # include .stan files of certain response distributions
  # Args:
  #   family: the model family
  # Returns:
  #   a list of character strings
  stopifnot(is(family, "family"))
  family <- family$family
  out <- list()
  if (family == "categorical") {
    out$data <- "  int<lower=2> ncat;  // number of categories \n" 
    out$tdataD <- "  vector[1] zero; \n"
    out$tdataC <- "  zero[1] = 0; \n"
  } else if (family == "zero_inflated_poisson") {
    out$fun <- "  #include 'fun_zero_inflated_poisson.stan' \n"
  } else if (family == "zero_inflated_negbinomial") {
    out$fun <- "  #include 'fun_zero_inflated_negbinomial.stan' \n"
  } else if (family == "zero_inflated_binomial") {
    out$fun <- "  #include 'fun_zero_inflated_binomial.stan' \n"
  } else if (family == "zero_inflated_beta") {
    out$fun <- "  #include 'fun_zero_inflated_beta.stan' \n"
  } else if (family == "hurdle_poisson") {
    out$fun <- "  #include 'fun_hurdle_poisson.stan' \n"
  } else if (family == "hurdle_negbinomial") {
    out$fun <- "  #include 'fun_hurdle_negbinomial.stan' \n"
  } else if (family == "hurdle_gamma") {
    out$fun <- "  #include 'fun_hurdle_gamma.stan' \n"
  } else if (family == "hurdle_lognormal") {
    out$fun <- "  #include 'fun_hurdle_lognormal.stan' \n"
  } else if (family == "exgaussian") {
    out$fun <- "  #include 'fun_exgaussian.stan' \n"
  } else if (family == "inverse.gaussian") {
    out$fun <- "  #include 'fun_inv_gaussian.stan' \n"
    out$tdataD <- "  #include 'tdataD_inv_gaussian.stan' \n"
    out$tdataC <- "  #include 'tdataC_inv_gaussian.stan' \n"
  } else if (family == "von_mises") {
    out$fun <- paste0("  #include 'fun_tan_half.stan' \n",
                      "  #include 'fun_von_mises.stan' \n")
  } else if (family == "wiener") {
    out$fun <- "  #include 'fun_wiener_diffusion.stan' \n"
    out$tdataD <- "  real min_Y; \n"
    out$tdataC <- "  min_Y = min(Y); \n"
  }
  out
}

stan_cens <- function(cens, family = gaussian()) {
  out <- list()
  if (cens) {
    stopifnot(is(family, "family"))
    out$data <- paste0(
      "  int<lower=-1,upper=2> cens[N];  // indicates censoring \n",
      if (isTRUE(attr(cens, "interval"))) {
        paste0(ifelse(use_int(family), " int rcens[N];", "  vector[N] rcens;"),
               "  // right censor points for interval censoring \n")
      })
  }
  out
}

stan_disp <- function(effects, family = gaussian()) {
  # stan code for models with addition argument 'disp'
  # Args:
  #   disp: logical; are dispersion factors specified?
  #   family: the model family
  stopifnot(is(family, "family"))
  out <- list()
  if (is(effects$disp, "formula")) {
    par <- if (has_sigma(family)) "sigma"
           else if (has_shape(family)) "shape"
           else stop("invalid family for addition argument 'disp'")
    if (!is.null(effects[[par]])) {
      stop2("Specifying 'disp' is not allowed when predicting '", par, "'.")
    }
    out$data <- "  vector<lower=0>[N] disp;  // dispersion factors \n"
    out$modelD <- paste0("  vector[N] disp_", par, "; \n")
    out$modelC1 <- paste0("  disp_", par, " = ", par, " * disp; \n")
  }
  out
}

stan_monotonic <- function(x) {
  # add the monotonic function to Stan's functions block
  if (grepl("[^[:alnum:]]monotonic\\(", collapse(x))) {
    "  #include fun_monotonic.stan \n"
  } else NULL
}

stan_misc_functions <- function(family = gaussian(), kronecker = FALSE) {
  # stan code for user defined functions
  # Args:
  #   family: the model family
  #   kronecker: logical; is the kronecker product needed?
  # Returns:
  #   a string containing defined functions in stan code
  stopifnot(is(family, "family"))
  out <- NULL
  if (family$link == "cauchit") {
    out <- paste0(out, "  #include 'fun_cauchit.stan' \n")
  } else if (family$link == "cloglog") {
    out <- paste0(out, "  #include 'fun_cloglog.stan' \n")
  }
  if (kronecker) {
    out <- paste0(out, "  #include 'fun_as_matrix.stan' \n",
                  "  #include 'fun_kronecker.stan' \n")
  }
  out
}

stan_prior <- function(class, coef = "", group = "", nlpar = "", suffix = "", 
                       prior = brmsprior(), wsp = 2, matrix = FALSE) {
  # Define priors for parameters in Stan language
  # Args:
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   nlpar: the name of a non-linear parameter
  #   suffix: a suffix to put at the parameter class
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   matrix: logical; corresponds the class to a parameter matrix?
  #   wsp: an integer >= 0 defining the number of spaces 
  #        in front of the output string
  # Returns:
  #   A character strings in stan language that defines priors 
  #   for a given class of parameters. If a parameter has has 
  #   no corresponding prior in prior, an empty string is returned.
  # only consider user defined priors related to this class and group
  wsp <- collapse(rep(" ", wsp))
  prior_only <- isTRUE(attr(prior, "prior_only"))
  hs_df <- attr(prior, "hs_df")
  keep <- prior$class == class & 
          (prior$coef %in% coef | !nzchar(prior$coef)) &
          (prior$group == group | !nzchar(prior$group)) & 
          (prior$nlpar %in% nlpar | !nzchar(prior$nlpar))
  prior <- prior[keep, ]
  if (!nchar(class) && nrow(prior)) {
    # increment_log_prob statements are directly put into the Stan code
    return(collapse(wsp, prior$prior, "; \n"))
  } 
  
  unlpar <- unique(nlpar)
  if (length(unlpar) > 1L) {
    # can only happen for SD parameters of the same ID
    base_prior <- rep(NA, length(unlpar))
    for (i in seq_along(unlpar)) {
      nlpar_prior <- prior[prior$nlpar %in% c("", unlpar[i]), ]
      base_prior[i] <- stan_base_prior(nlpar_prior)
    }
    if (length(unique(base_prior)) > 1L) {
      # define prior for single coefficients manually
      # as there is not single base_prior anymore
      take <- match(prior[nzchar(prior$coef), "nlpar"], unlpar)
      prior[nzchar(prior$coef), "prior"] <- base_prior[take]
    }
    base_prior <- base_prior[1]
  } else {
    base_prior <- stan_base_prior(prior)
  }
  
  individual_prior <- function(i, max_index) {
    # individual priors for each parameter of a class
    if (max_index > 1L || matrix) {
      index <- paste0("[",i,"]")      
    } else {
      index <- ""
    }
    if (length(nlpar) > 1L) {
      prior <- prior[prior$nlpar == nlpar[i], ]
    }
    uc_prior <- prior$prior[match(coef[i], prior$coef)]
    if (!is.na(uc_prior) & nchar(uc_prior)) { 
      # user defined prior for this parameter
      coef_prior <- uc_prior
    } else { # base prior for this parameter
      coef_prior <- base_prior 
    }  
    if (nchar(coef_prior) > 0) {  # implies a proper prior
      out <- paste0(wsp, class, index, " ~ ", coef_prior, "; \n")
    } else {
      out <- ""  # implies an improper flat prior
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(class, suffix)
  if (any(with(prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(seq_along(coef), individual_prior, 
                  max_index = length(coef))
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- paste0("to_vector(", class, ")")
    }
    out <- paste0(wsp, class, " ~ ", base_prior, "; \n")
  } else {
    out <- ""
  }
  if (class == "b" && !is.null(hs_df)) {
    # add horseshoe shrinkage priors
    hs_shrinkage_priors <- paste0(
      "  hs_local ~ student_t(", attr(prior, "hs_df"), ", 0, 1); \n",
      "  hs_global ~ cauchy(0, ", attr(prior, "hs_scale_global"), "); \n")
    out <- c(hs_shrinkage_priors, out)
  }
  out <- collapse(out)
  if (prior_only && nchar(class) && !nchar(out)) {
    stop2("Sampling from priors is not possible because ", 
          "not all parameters have proper priors. \n",
          "Error occured for class '", class, "'.")
  }
  out
}

stan_base_prior <- function(prior) {
  # get base (highest level) prior of all priors
  # Args:
  #   prior: a prior.frame
  stopifnot(length(unique(prior$class)) <= 1L) 
  igroup <- which(with(prior, !nchar(coef) & nchar(group) & nchar(prior)))
  inlpar <- which(with(prior, !nchar(coef) & nchar(nlpar) & nchar(prior)))
  iclass <- which(with(prior, !nchar(coef) & !nchar(group) & nchar(prior)))
  if (length(igroup)) {  
    # if there is a global prior for this group
    base_prior <- prior[igroup, "prior"]
  } else if (length(inlpar)) {
    # if there is a global prior for this non-linear parameter
    base_prior <- prior[inlpar, "prior"]
  } else if (length(iclass)) {  
    # if there is a global prior for this class
    base_prior <- prior[iclass, "prior"]
  } else {  
    # no proper prior for this class
    base_prior <- ""
  }
  stopifnot(length(base_prior) == 1L)
  base_prior
}

stan_rngprior <- function(sample_prior, prior, par_declars = "",
                          family = gaussian(), hs_df = NULL) {
  # stan code to sample from priors seperately
  #
  # Args:
  #   sample_prior: take samples from priors?
  #   prior: character string taken from stan_prior
  #   par_declars: the parameters block of the Stan code
  #                requied to extract boundaries
  #   family: the model family
  #   hs_df: hs_df degrees of freedom
  #
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  stopifnot(is(family, "family"))
  out <- list()
  if (sample_prior) {
    prior <- gsub(" ", "", paste0("\n", prior))
    pars <- gsub("\\\n|to_vector\\(|\\)", "", get_matches("\\\n[^~]+", prior))
    take <- !grepl("^(z|zs)_|^increment_log_prob\\(", pars)
    pars <- rename(pars[take], symbols = c("^L_", "^Lrescor"), 
                   subs = c("cor_", "rescor"), fixed = FALSE)
    dis <- gsub("~", "", regmatches(prior, gregexpr("~[^\\(]+", prior))[[1]])[take]
    args <- regmatches(prior, gregexpr("\\([^;~]+\\);", prior))[[1]][take]
    type <- rep("real", length(pars))
    
    # rename parameters containing indices
    has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
    pars[has_ind] <- ulapply(pars[has_ind], function(par) {
      ind <- regmatches(par, gregexpr("\\[[[:digit:]]+\\]", par))
      ind <- as.numeric(substr(ind, 2, nchar(ind) - 1))
      gsub("\\[[[:digit:]]+\\]", paste0("_", ind), par)
    })
    
    # special treatment of lkj_corr_cholesky priors
    args <- ifelse(grepl("corr_cholesky$", dis), 
                   paste0("(2,", substr(args, 2, nchar(args)-1), "[1, 2];"), 
                   args)
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
      out$par <- paste0("  // parameters to store prior samples \n",
                        collapse("  real", bounds[has_bounds], 
                                 " prior_", pars[has_bounds], "; \n"))
      out$model <- paste0("  // additionally draw samples from priors \n",
                          collapse("  prior_", pars[has_bounds] ," ~ ",
                                   dis[has_bounds], args[has_bounds]," \n"))
    }
    no_bounds <- !has_bounds
    if (any(no_bounds)) {  
      # unbounded parameters can be sampled in the generatated quantities block
      if (!is.null(hs_df)) {
        args[grepl("^b(m|p|_|$)", pars)] <- "(0, prior_hs_local * prior_hs_global);" 
      }
      out$genD <- collapse(
        "  ", types[no_bounds], " prior_", pars[no_bounds], "; \n")
      out$genC <- paste0(
        "  // additionally draw samples from priors \n",
        collapse("  prior_", pars[no_bounds], " = ",
                 dis[no_bounds], "_rng", args[no_bounds], " \n"))
    }
    # compute priors for the actual population-level intercepts
    is_temp_intercept <- grepl("^temp.*_Intercept", pars)
    if (any(is_temp_intercept)) {
      temp_intercepts <- pars[is_temp_intercept]
      p <- gsub("^temp|_Intercept$", "", temp_intercepts)
      intercepts <- paste0("b", p, "_Intercept")
      use_plus <- family$family %in% c("cumulative", "sratio")
      sub_X_means <- paste0(ifelse(use_plus, " + ", " - "), 
                            "dot_product(means_X", p, ", b", p, ")")
      out$genD <- paste0(out$genD, 
        collapse("  real prior_", intercepts, "; \n"))
      out$genC <- paste0(out$genC, 
        collapse("  prior_", intercepts, " = ",
                 "prior_", temp_intercepts, sub_X_means, "; \n"))
    }
  }
  out
}

stan_link <- function(link) {
  # find the link in Stan language
  # Args:
  #   link: the link function
  switch(link, identity = "", log = "log", inverse = "inv",
         sqrt = "sqrt", "1/mu^2" = "inv_square", logit = "logit", 
         probit = "inv_Phi", probit_approx = "inv_Phi", 
         cloglog = "cloglog", cauchit = "cauchit",
         tan_half = "tan_half")
}

stan_ilink <- function(link) {
  # find the inverse link in Stan language
  # Args:
  #   link: the link function
  switch(link, identity = "", log = "exp", inverse = "inv", 
         sqrt = "square", "1/mu^2" = "inv_sqrt", logit = "inv_logit", 
         probit = "Phi", probit_approx = "Phi_approx", 
         cloglog = "inv_cloglog", cauchit = "inv_cauchit",
         tan_half = "inv_tan_half")
}

stan_ll_adj <- function(effects, adds = c("weights", "cens", "trunc")) {
  # checks if certain 'adds' are present so that the LL has to be adjusted
  # Args:
  #   effects: output of extract_effects
  #   adds: vector of addition argument names
  stopifnot(all(adds %in% c("weights", "cens", "trunc")))
  any(ulapply(effects[adds], is.formula))
}

stan_has_built_in_fun <- function(family, link) {
  # indicates if a family-link combination has a build in 
  # function in Stan (such as binomial_logit)
  (family %in% c("binomial", "bernoulli", "cumulative", "categorical")
   && link == "logit" || is.count(family) && link == "log")
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

stan_apn <- function(ap, auxpars, reqn = NULL) {
  # helper function to add "[n]" to auxiliary parameter strings
  paste0(ap, if ((is.null(reqn) || reqn) && ap %in% auxpars) "[n]")
}
