stan_fixef <- function(fixef, intercepts = "Intercept", 
                       family = gaussian(), prior = prior_frame(), 
                       sparse = FALSE, threshold = "flexible") {
  # Stan code for fixed effects
  #
  # Args:
  #   fixef: names of the fixed effects
  #   csef: names of the category specific effects
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   nint: number of fixed effects intercepts
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   a list containing Stan code related to fixed effects
  nint <- length(intercepts)
  out <- list()
  if (length(fixef)) {
    centered <- ifelse(nint > 0, "centered", "")
    out$data <- paste0(out$data, 
      "  int<lower=1> K;  // number of population-level effects \n", 
      "  matrix[N, K] X;  // ", centered, " population-level design matrix \n")
    if (sparse) {
      out$tdataD <- "  #include tdata_def_sparse_X.stan \n"
      out$tdataC <- "  #include tdata_calc_sparse_X.stan \n"
    }
    bound <- with(prior, bound[class == "b" & coef == ""])
    out$par <- paste0(out$par,
      "  vector", bound, "[K] b;  // population-level effects \n") 
    fixef_prior <- stan_prior(class = "b", coef = fixef, prior = prior)
    out$prior <- paste0(out$prior, fixef_prior)
  }
  if (nint) {
    if (length(fixef)) {
      # X_means is not defined in intercept only models
      def_X_means <- paste0("  vector[K] X_means", 
                            if (nint > 1L) "[nint]", 
                            ";  // column means of X before centering \n")
      sub_X_means <- paste0(" - dot_product(X_means", 
                            if (nint > 1L) "[i]", ", b)")
    } else {
      def_X_means <- sub_X_means <- ""
    }
    if (is.ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$data <- paste0(out$data, def_X_means)
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      out$genC <- paste0("  b_Intercept <- temp_Intercept", sub_X_means, "; \n")
    } else {
      if (nint == 1L) {
        out$data <- paste0(out$data, def_X_means)
        out$par <- paste0(out$par, 
          "  real temp_Intercept;  // temporary Intercept \n")
        out$genD <- "  real b_Intercept;  // population-level intercept \n"
        out$genC <- paste0("  b_Intercept <- temp_Intercept", sub_X_means, "; \n")
      } else if (nint > 1L) {
        out$data <- paste0(out$data,
          "  int nint;  // number of population-level intercepts \n",
          "  int Jint[N];  // assigns intercepts to observations \n", 
          def_X_means)
        out$par <- paste0(out$par, 
          "  vector[nint] temp_Intercept;  // temporary Intercepts \n")
        out$genD <- "  vector[nint] b_Intercept;  // population-level intercepts \n"
        out$genC <- paste0(
          "  for (i in 1:nint) { \n",
          "    b_Intercept[i] <- temp_Intercept[i]", sub_X_means, "; \n",
          "  } \n")
      }
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    suffix <- ifelse(threshold == "equidistant", "1", "")
    int_prior <- stan_prior("temp_Intercept", coef = intercepts, 
                            prior = prior, suffix = suffix)
    out$prior <- paste0(out$prior, int_prior)
  }
  out
}

stan_ranef <- function(i, ranef, prior = prior_frame(), 
                       names_cov_ranef = NULL, nlpar = NULL) {
  # Random effects in Stan 
  # 
  # Args:
  #   i: the index of the grouping factor
  #   ranef: a named list returned by gather_ranef
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   names_cov_ranef: names of the grouping factors 
  #                    for which custom covariance matrices are specified.
  #   par: an optional character string to add to the variable names
  #        (used for non-linear models)
  #
  # Returns:
  #   A vector of strings containing the random effects in stan language
  r <- ranef[[i]]
  g <- attr(ranef[[i]], "group")
  cor <- attr(ranef[[i]], "cor")
  ccov <- g %in% names_cov_ranef
  if (!is.null(nlpar)) {
    i <- paste0(nlpar, "_", i)
  }
  out <- list()
  out$data <- paste0(
    "  // data for group-specific effects of ", g, " \n",
    "  int<lower=1> J_", i, "[N]; \n",
    "  int<lower=1> N_", i, "; \n",
    "  int<lower=1> K_", i, "; \n",
    if (ccov) paste0(
      "  // cholesky factor of known covariance matrix \n",
      "  matrix[N_", i, ", N_", i,"] Lcov_", i,"; \n"))
  
  out$prior <- stan_prior(class = "sd", group = g, coef = r, nlpar = nlpar,
                          suffix = paste0("_", i), prior = prior)
  if (length(r) == 1L) {  # only one random effect
    out$data <- paste0(out$data, 
      "  vector[N] Z_", i, "; \n")
    out$par <- paste0(
      "  real<lower=0> sd_", i, ";",
      "  // group-specific standard deviation \n",
      "  vector[N_", i, "] z_", i, ";",
      "  // unscaled group-specific effects \n")
    out$prior <- paste0(out$prior,"  z_", i, " ~ normal(0, 1); \n")
    out$transD <- paste0(
      "  // group-specific effects \n",
      "  vector[N_", i, "] r_", i, "; \n")
    out$transC <- paste0("  r_", i,  " <- sd_", i, " * (", 
                         if (ccov) paste0("Lcov_", i, " * "), "z_", i, ");\n")
  } else if (length(r) > 1L) {
    j <- seq_along(r)
    out$data <- paste0(out$data, 
      collapse("  vector[N] Z_", i, "_", j, ";  \n"))
    out$par <- paste0(
      "  vector<lower=0>[K_", i, "] sd_", i, ";",
      "  // group-specific standard deviations \n")
    if (cor) {  
      # multiple correlated random effects
      out$data <- paste0(out$data,  
        "  int<lower=1> NC_", i, "; \n")
      out$par <- paste0(out$par,
        "  matrix[K_", i, ", N_", i, "] z_", i, ";",
        "  // unscaled group-specific effects \n",    
        "  // cholesky factor of correlation matrix \n",
        "  cholesky_factor_corr[K_", i, "] L_", i, "; \n")
      out$prior <- paste0(out$prior, 
        stan_prior(class = "L", group = g, nlpar = nlpar,
                   suffix = paste0("_", i), prior = prior),
        "  to_vector(z_", i, ") ~ normal(0, 1); \n")
      out$transD <- paste0(
        "  // group-specific effects \n",
        "  matrix[N_", i, ", K_", i, "] r_", i, "; \n",
        collapse("  vector[N_", i, "] r_", i, "_", j, "; \n"))
      if (ccov) {  # customized covariance matrix supplied
        out$transC <- paste0("  r_", i," <- as_matrix(kronecker(Lcov_", i, ",", 
          " diag_pre_multiply(sd_", i,", L_", i,")) *",
          " to_vector(z_", i, "), N_", i, ", K_", i, "); \n")
      } else { 
        out$transC <- paste0("  r_", i, " <- ", 
          "(diag_pre_multiply(sd_", i, ", L_", i,") * z_", i, ")'; \n")
      }
      out$transC <- paste0(out$transC, 
        collapse("  r_", i, "_", j, " <- r_", i, "[, ",j,"];  \n"))
      # return correlations above the diagonal only
      cors_genC <- ulapply(2:length(r), function(k) 
        lapply(1:(k - 1), function(j) paste0(
          "  cor_", i, "[", (k - 1) * (k - 2) / 2 + j, 
          "] <- Cor_", i, "[", j, ",", k, "]; \n")))
      out$genD <- paste0(
        "  corr_matrix[K_", i, "] Cor_", i, "; \n",
        "  vector<lower=-1,upper=1>[NC_", i, "] cor_", i, "; \n")
      out$genC <- paste0(
        "  // take only relevant parts of correlation matrix \n",
        "  Cor_", i, " <- multiply_lower_tri_self_transpose(L_", i, "); \n",
        collapse(cors_genC)) 
    } else {
      # multiple uncorrelated random effects
      out$par <- paste0(out$par,
        "  vector[N_", i, "] z_", i, "[K_", i, "];",
        "  // unscaled group-specific effects \n")
      out$prior <- paste0(out$prior, collapse(
        "  z_", i, "[", j, "] ~ normal(0, 1); \n"))
      out$transD <- paste0("  // group-specific effects \n", 
        collapse("  vector[N_", i, "] r_", i, "_", j, "; \n"))
      out$transC <- collapse(
        "  r_", i, "_", j, " <- sd_", i, "[", j, "] * (", 
        if (ccov) paste0("Lcov_", i, " * "), "z_", i, "[", j, "]); \n")
      out$genD <- paste0(
        "  matrix[N_", i, ", K_", i, "] r_", i, "; \n")
      out$genC <- collapse(
        "  r_", i, "[, ", j, "] <- r_", i, "_", j, "; \n")
    }
  }
  out
}

stan_monef <- function(monef, prior = prior_frame()) {
  # Stan code for monotonous effects
  # Args:
  #   csef: names of the monotonous effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  out <- list()
  if (length(monef)) {
    I <- seq_along(monef)
    out$fun <- "  #include fun_monotonous.stan \n"
    out$data <- paste0(
      "  int<lower=1> Km; \n",
      "  int Xm[N, Km]; \n",
      "  int<lower=2> Jm[Km]; \n",
      collapse("  vector[Jm[", I, "]] con_simplex_", I, "; \n"))
    bound <- with(prior, bound[class == "b" & coef == ""])
    out$par <- paste0(
      "  vector", bound, "[Km] bm; \n",
      collapse("  simplex[Jm[", I, "]] simplex_", I, "; \n")) 
    out$prior <- paste0(
      stan_prior(class = "b", coef = monef, prior = prior, suffix = "m"),
      collapse("  simplex_", I, " ~ dirichlet(con_simplex_", I, "); \n"))
  }
  out
}

stan_csef <- function(csef, prior = prior_frame()) {
  # Stan code for category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  out <- list()
  if (length(csef)) {
    out$data <- paste0(
      "  int<lower=1> Kp;  // number of category specific effects \n",
      "  matrix[N, Kp] Xp;  // CSE design matrix \n")
    bound <- with(prior, bound[class == "b" & coef == ""])
    out$par <- paste0(
      "  matrix", bound, "[Kp, ncat - 1] bp;  // category specific effects \n")
    out$prior <- stan_prior(class = "b", coef = csef, prior = prior, 
                             suffix = "p", matrix = TRUE)
  }
  out
} 

stan_llh <- function(family, se = FALSE, weights = FALSE, trials = FALSE, 
                     cens = FALSE, disp = FALSE, trunc = .trunc(), 
                     autocor = cor_arma(), cse = FALSE, nresp = 1L) {
  # Likelihoods in stan language
  #
  # Args:
  #   family: the model family
  #   se: logical; user defined SEs present?
  #   weights: logical; weights present?
  #   trials: logical; number of bernoulli trials given per observation?
  #   cens: logical; censored data?
  #   trunc: list containing lower and upper truncation boundaries
  #   autocor: autocorrelation structure; an object of classe cor_arma
  #   cse: a flag indicating if category specific effects are present
  #        (for ordinal models only)
  #   nresp: number of response variables
  #
  # Returns:
  #   a string containing the likelihood of the model in stan language
  stopifnot(is(family, "family"))
  link <- family$link
  type <- family$type
  family <- family$family
  is_categorical <- is.categorical(family)
  is_ordinal <- is.ordinal(family)
  is_hurdle <- is.hurdle(family)
  is_zero_inflated <- is.zero_inflated(family)
  is_forked <- is.forked(family)
  is_multi <- is.linear(family) && nresp > 1L
  is_trunc <- trunc$lb > -Inf || trunc$ub < Inf
  if (is_multi) {
    # prepare for use of a multivariate likelihood
    family <- paste0("multi_", family)
  } else if (use_cov(autocor) && (get_ar(autocor) || get_ma(autocor))) {
    # ARMA effects have a special formulation
    # if fitted using a covariance matrix for residuals
    family <- paste0(family, "_cov")
    if (weights || cens || is_trunc) {
      stop("Invalid addition arguments", call. = FALSE)
    }
  } else if (is.lognormal(family, link = link)) {
    # prepare for use of lognormal likelihood
    family <- "lognormal"
    link <- "identity"
  }
  
  simplify <- !is_trunc && !cens && has_build_in_fun(family, link)
  n <- ifelse(cens || weights || is_trunc || is_categorical ||
              is_ordinal || is_hurdle || is_zero_inflated, "[n]", "")
  ns <- ifelse((se || trials || disp) && (cens || weights || is_trunc) 
               || (trials && is_zero_inflated), "[n]", "")
  disp <- ifelse(disp, "disp_", "")
  sigma <- paste0(ifelse(se, "se", paste0(disp, "sigma")), ns)
  shape <- paste0(disp, "shape", ns)
  ordinal_args <- paste0("eta[n], ", if (cse) "etap[n], ", 
                         "temp_Intercept")
  
  # use inverse link in likelihood statement only 
  # if it does not prevent vectorization 
  ilink <- ifelse(n == "[n]" && !simplify, stan_ilink(link), "")
  eta <- paste0("eta", if (!is.null(type)) paste0("_", type))
  if (n == "[n]") {
    if (is_hurdle || is_zero_inflated) {
      eta <- paste0(ilink,"(eta[n]), ", ilink,"(eta[n + N_trait])")
    } else {
      eta <- paste0(eta, "[n]")
      fl <- ifelse(family %in% c("gamma", "exponential"), 
                   paste0(family,"_",link), family)
      eta <- switch(fl, paste0(ilink,"(",eta,")"),
                    gamma_log = paste0(shape, " * exp(-",eta,")"),
                    gamma_inverse = paste0(shape, " * ", eta),
                    gamma_identity =  paste0(shape, " / ", eta),
                    exponential_log = paste0("exp(-",eta,")"),
                    exponential_inverse = eta,
                    exponential_identity = paste0("inv(",eta,")"),
                    weibull = paste0(ilink, "(",eta, " / ", shape,")"))
    }
  }
  
  if (simplify) { 
    llh_pre <- switch(family,
      poisson = c("poisson_log", eta), 
      negbinomial = c("neg_binomial_2_log", paste0(eta,", ",shape)),
      geometric = c("neg_binomial_2_log", paste0(eta,", 1")),
      cumulative = c("ordered_logistic", "eta[n], temp_Intercept"),
      categorical = c("categorical_logit", "append_row(zero, eta[J_trait[n]])"),
      binomial = c("binomial_logit", paste0("trials",ns,", ",eta)), 
      bernoulli = c("bernoulli_logit", eta))
  } else {
    llh_pre <- switch(family,
      gaussian = c("normal", paste0(eta,", ",sigma)),
      gaussian_cov = c("normal_cov", paste0(eta,", se2, N_tg, ", 
                       "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      student = c("student_t",  paste0("nu, ",eta,", ",sigma)),
      student_cov = c("student_t_cov", paste0("nu, ",eta,", se2, N_tg, ", 
                      "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      cauchy = c("cauchy", paste0(eta,", ", sigma)),
      cauchy_cov = c("student_t_cov", paste0("1, ",eta,", se2, N_tg, ", 
                     "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      lognormal = c("lognormal", paste0(eta,", ",sigma)),
      multi_gaussian = c("multi_normal_cholesky", paste0("Eta",n,", LSigma")),
      multi_student = c("multi_student_t", paste0("nu, Eta",n,", Sigma")),
      multi_cauchy = c("multi_student_t", paste0("1.0, Eta",n,", Sigma")),
      poisson = c("poisson", eta),
      negbinomial = c("neg_binomial_2", paste0(eta,", ",shape)),
      geometric = c("neg_binomial_2", paste0(eta,", 1")),
      binomial = c("binomial", paste0("trials",ns,", ",eta)),
      bernoulli = c("bernoulli", eta),
      gamma = c("gamma", paste0(shape,", ",eta)), 
      exponential = c("exponential", eta),
      weibull = c("weibull", paste0(shape,", ",eta)), 
      inverse.gaussian = c(paste0("inv_gaussian", if (!nchar(n)) "_vector"), 
                           paste0(eta, ", shape, log_Y",n,", sqrt_Y",n)),
      beta = c("beta", paste0(eta, " * phi, (1 - ", eta, ") * phi")),
      cumulative = c("cumulative", ordinal_args),
      sratio = c("sratio", ordinal_args),
      cratio = c("cratio", ordinal_args),
      acat = c("acat", ordinal_args),
      hurdle_poisson = c("hurdle_poisson", "eta[n], eta[n + N_trait]"),
      hurdle_negbinomial = c("hurdle_neg_binomial_2", 
                             "eta[n], eta[n + N_trait], shape"),
      hurdle_gamma = c("hurdle_gamma", "shape, eta[n], eta[n + N_trait]"),
      zero_inflated_poisson = c("zero_inflated_poisson", 
                                "eta[n], eta[n + N_trait]"),
      zero_inflated_negbinomial = c("zero_inflated_neg_binomial_2", 
                                    "eta[n], eta[n + N_trait], shape"),
      zero_inflated_binomial = c("zero_inflated_binomial", 
         paste0("trials",ns,", eta[n], eta[n + N_trait]")),
      zero_inflated_beta = c("zero_inflated_beta", 
                             "eta[n], eta[n + N_trait], phi"))
  }
  
  # write likelihood code
  type <- c("cens", "weights")[match(TRUE, c(cens, weights))]
  if (is.na(type)) type <- "general"
  # prepare for possible truncation
  code_trunc <- ""
  if (is_trunc) {
    if (type %in% c("cens", "weights")) {
      stop("truncation is not yet possible in censored or weighted models",
           call. = FALSE)
    } else {
      lb <- ifelse(trunc$lb > -Inf, "lb", "")
      ub <- ifelse(trunc$ub < Inf, "ub", "")
      code_trunc <- paste0(" T[",lb,", ",ub,"]")
    }
  }
  add_weights <- ifelse(weights, "weights[n] * ", "")
  llh <- switch(type, 
    cens = paste0("  // special treatment of censored data \n",
      "    if (cens[n] == 0) ", 
      ifelse(!weights, paste0("Y[n] ~ ", llh_pre[1],"(", llh_pre[2],"); \n"),
             paste0("increment_log_prob(", add_weights, 
                    llh_pre[1], "_log(Y[n], ", llh_pre[2],")); \n")),
      "    else { \n",         
      "      if (cens[n] == 1) increment_log_prob(", add_weights, 
      llh_pre[1], "_ccdf_log(Y[n], ", llh_pre[2],")); \n",
      "      else increment_log_prob(", add_weights, 
      llh_pre[1], "_cdf_log(Y[n], ", llh_pre[2],")); \n",
      "    } \n"),
    weights = paste0("  lp_pre[n] <- ", llh_pre[1], "_log(Y[n], ",
                     llh_pre[2],"); \n"),
    general = paste0("  Y", n, " ~ ", llh_pre[1],"(", llh_pre[2],")", 
                     code_trunc, "; \n")) 
  # loop over likelihood if it cannot be vectorized
  trait <- ifelse(is_multi || is_forked || is_categorical, "_trait", "")
  if (weights || cens || is_trunc || is_ordinal || is_categorical || 
      is_hurdle || is_zero_inflated) {
    llh <- paste0("  for (n in 1:N", trait, ") { \n  ", llh, "  } \n")
  }
  llh
}

stan_eta <- function(family, fixef, ranef = list(), csef = NULL, 
                     monef = NULL, nint = 1, autocor = cor_arma(),  
                     sparse = FALSE, add = FALSE, disp = FALSE, 
                     offset = FALSE, is_multi = FALSE) {
  # linear predictor in Stan
  #
  # Args:
  #   family: the model family
  #   fixef: names of the fixed effects
  #   ranef: a named list returned by gather_ranef
  #   csef: names of the category specific effects
  #   nint: number of fixed effects intercepts
  #   autocor: autocorrelation structure
  #   sparse: is the fixed effects design matrix sparse?
  #   add: is the model weighted, censored, or truncated?
  #   disp: is the 'disp' addition argument specified?
  #   offset: is an offset defined?
  #   is_multi: is the model multivariate?
  # 
  # Return:
  #   the linear predictor in stan language
  stopifnot(is(family, "family"))
  link <- family$link
  family <- family$family
  # initialize eta
  eta <- list()
  eta$transD <- paste0(
    "  vector[N] eta;  // linear predictor \n", 
    if (length(csef)) 
      paste0("  matrix[N, ncat - 1] etap;",
             "  // linear predictor for category specific effects \n"),
    if (is_multi) 
      paste0("  vector[K_trait] Eta[N_trait];",
             "  // multivariate linear predictor matrix \n"))
  eta_obj <- ifelse(is_multi, "Eta[m, k]", "eta[n]")
  wsp <- ifelse(is_multi, "  ", "")
  
  # transform eta before it is passed to the likelihood
  eta$transform <- stan_eta_transform(family, link, add = add)
  eta_ilink <- rep("", 2)
  if (eta$transform || (get_ar(autocor) && !use_cov(autocor))) {
    eta_ilink <- stan_eta_ilink(family, link, disp = disp)
    if (get_ar(autocor)) {
      eta_ar <- ifelse(!use_cov(autocor), " + head(E[n], Kar) * ar", "")
      eta$transC3 <- paste0("    ", wsp, eta_obj," <- ", eta_ilink[1], 
                            eta_obj, eta_ar, eta_ilink[2], "; \n")
      eta_ilink <- rep("", 2)  # don't apply the link function twice
    }
  }
  
  # define fixed, random, and autocorrelation effects
  eta_int <- if (nint > 1L) " + temp_Intercept[Jint[n]]"
  eta_ranef <- stan_eta_ranef(ranef)
  eta_monef <- stan_eta_monef(monef)
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor), 
                   " + head(E[n], Kma) * ma", "")
  add2eta <- any(nchar(c(eta_int, eta_monef, eta_ma, eta_ranef, eta_ilink[1])))
  if (add2eta || is_multi) {
    eta$transC2 <- paste0(
      "    ", wsp, eta_obj," <- ", eta_ilink[1], "eta[n]", 
      eta_int, eta_monef, eta_ranef, eta_ma, eta_ilink[2],"; \n")
  }
  eta_fixef <- stan_eta_fixef(fixef, sparse = sparse)
  eta_cse <- if (length(csef)) "  etap <- Xp * bp; \n"
  eta$transC1 <- paste0(
    "  // compute linear predictor \n",
    "  eta <- ", eta_fixef, 
    if (nint == 1L && !is.ordinal(family)) " + temp_Intercept",
    if (offset) " + offset",
    if (get_arr(autocor)) " + Yarr * arr", 
    "; \n", eta_cse)
  eta
}

stan_nonlinear <- function(effects, data, family = gaussian(), 
                           add = FALSE, cov_ranef = NULL, 
                           prior = prior_frame(), disp = FALSE) {
  # prepare Stan code for non-linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   add: Is the model weighted, censored, or truncated?
  #   cov_ranef: a list of user-defined covariance matrices
  #   prior: a prior_frame object
  #   disp: is the 'disp' addition argument specified?
  out <- list()
  if (length(effects$nonlinear)) {
    for (i in seq_along(effects$nonlinear)) {
      nlp <- names(effects$nonlinear)[i]
      eta <- paste0("eta_", nlp)
      out$transD <- paste0(out$transD, "  vector[N] ", eta, "; \n")
      out$data <- paste0(out$data, 
        "  // data for non-linear effects of ", nlp, "\n")
      out$par <- paste0(out$par, 
        "  // non-linear effects of ", nlp, "\n")
      # include fixed effects
      fixef <- colnames(get_model_matrix(effects$nonlinear[[i]]$fixed, data))
      if (length(fixef)) {
        out$data <- paste0(out$data, 
          "  int<lower=1> K_", nlp, "; \n", 
          "  matrix[N, K_", nlp, "] X_", nlp, "; \n")
        bound <- with(prior, bound[class == "b" & coef == "" & nlpar == nlp])
        out$par <- paste0(out$par,
         "  vector", bound, "[K_", nlp, "] b_", nlp, "; \n")
        out$transC1 <- paste0(out$transC1, 
          "  ", eta, " <- X_", nlp, " * b_", nlp, "; \n")  
        out$prior <- paste0(out$prior,
          stan_prior(class = "b", coef = fixef, nlpar = nlp, 
                     suffix = paste0("_", nlp), prior = prior))
      } else {
        out$transC1 <- paste0(out$transC1, 
          "  ", eta, " <- rep_vector(0, N); \n")  
      }
      # include random effects
      ranef <- gather_ranef(effects$nonlinear[[i]], data = data)
      if (length(ranef)) {
        text_ranef <- lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
                             names_cov_ranef = names(cov_ranef), 
                             prior = prior, nlpar = nlp)
        text_ranef <- collapse_lists(text_ranef)
        out$data <- paste0(out$data, text_ranef$data)
        out$prior <- paste0(out$prior, text_ranef$prior)
        out$par <- paste0(out$par, text_ranef$par)
        out$transD <- paste0(out$transD, text_ranef$transD)
        out$transC1 <- paste0(out$transC1, text_ranef$transC)
        out$transC2 <- paste0(out$transC2, 
          "    ", eta, "[n] <- ", eta, "[n]", 
          stan_eta_ranef(ranef, nlpar = nlp), "; \n") 
        out$genD <- paste0(out$genD, text_ranef$genD)
        out$genC <- paste0(out$genC, text_ranef$genC)
      }
      
    }
    # prepare non-linear model of eta 
    nlpars <- wsp(names(effects$nonlinear))
    new_nlpars <- paste0(" eta_", names(effects$nonlinear), "[n] ")
    # covariates in the nonlinear model
    covars <- wsp(setdiff(all.vars(effects$fixed[[3]]), 
                          names(effects$nonlinear)))
    if (length(covars)) {
      out$data <- paste0(out$data, 
        "  int<lower=1> KC;  // number of covariates \n",
        "  matrix[N, KC] C;  // covariate matrix \n")
      new_covars <- paste0(" C[n, ", seq_along(covars), "] ")
    } else new_covars <- NULL
    # add whitespaces to be able to replace parameters and covariates
    meta_sym <- c("+", "-", "*", "/", "^", ")", "(", ",")
    nlmodel <- gsub(" ", "", collapse(deparse(effects$fixed[[3]])))
    nlmodel <- wsp(rename(nlmodel, meta_sym, wsp(meta_sym))) 
    nlmodel <- rename(nlmodel, c(nlpars, covars, " ( ", " ) "), 
                      c(new_nlpars, new_covars, "(", ")"))
    # possibly transform eta in the transformed params block
    transform <- stan_eta_transform(family$family, family$link, add = add)
    if (transform) {
      eta_ilink <- stan_eta_ilink(family$family, family$link, disp = disp)
    } else eta_ilink <- rep("", 2)
    out$transD <- paste0(out$transD, "  vector[N] eta; \n")
    out$transC2 <- paste0(out$transC2, 
      "    // compute non-linear predictor \n",
      "    eta[n] <- ", eta_ilink[1], trimws(nlmodel), eta_ilink[2], "; \n")
  }
  out
}

stan_arma <- function(family, autocor, prior = prior_frame(),
                      has_se = FALSE, has_disp = FALSE, 
                      is_multi = FALSE, nonlinear = NULL) {
  # AR(R)MA autocorrelation in Stan
  # 
  # Args:
  #   family: the model family
  #   autocor: autocorrelation structure; object of class cor_arma
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   has_se: user defined standard errors present?
  #   is_multi: is the model multivariate?
  #   nonlinear: optional list of nonlinear formulas
  #
  # Returns:
  #   stan code for computing AR(R)MA effects
  stopifnot(is(family, "family"))
  is_linear <- is.linear(family)
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  Karr <- get_arr(autocor)
  out <- list()
  if (Kar || Kma) {
    if (!is_linear) {
      stop(paste("ARMA effects for family", family$family, 
                 "are not yet implemented"), call. = FALSE)
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
      err_msg <- "ARMA covariance matrices are not yet allowed"
      if (is_multi) {
        stop(paste(err_msg, "in multivariate models."), call. = FALSE)
      }
      if (has_disp) {
        stop(paste(err_msg, "when specifying 'disp'."), call. = FALSE)
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
                            "  res_cov_matrix <- cov_matrix_", cov_mat_fun, 
                            "(", cov_mat_args, ", sigma, max(nobs_tg)); \n")
      # defined selfmade functions for the functions block
      if (family$family == "gaussian") {
        out$fun <- paste0(out$fun, "  #include 'fun_normal_cov.stan' \n")
      } else { # family %in% c("student", "cauchy")
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
      if (has_se) {
        stop(paste(err_msg, "when specifying 'se'."), call. = FALSE)
      }
      if (length(nonlinear)) {
        stop(paste(err_msg, "for non-linear models."), call. = FALSE)
      }
      index <- ifelse(is_multi, "m, k", "n")
      wsp <- ifelse(is_multi, "      ", "    ")
      link <- c(identity = "", log = "log", inverse = "inv")[family$link]
      out$transD <- paste0("  matrix[N, Karma] E;  // ARMA design matrix \n",
                           "  vector[N] e;  // residuals \n") 
      out$transC1 <- "  E <- E_pre; \n" 
      out$transC2 <- paste0(
        wsp, "// calculation of ARMA effects \n",
        wsp, "e[n] <- ", link, "(Y[", index, "]) - eta[n]", "; \n",
        wsp, "for (i in 1:Karma) { \n", 
        wsp, "  if (n + 1 - i > 0 && n < N && tg[n + 1] == tg[n + 1 - i]) { \n",
        wsp, "     E[n + 1, i] <- e[n + 1 - i]; \n",
        wsp, "  } \n",
        wsp, "} \n")
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
  out
}

stan_multi <- function(family, response, prior = prior_frame()) {
  # some Stan code for multivariate models
  #
  # Args:
  #   family: model family
  #   response: names of the response variables
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # 
  # Returns: 
  #   list containing Stan code specific for multivariate models
  stopifnot(is(family, "family"))
  out <- list()
  nresp <- length(response)
  if (nresp > 1L) {
    if (is.linear(family)) {
      out$data <- "  #include 'data_multi.stan' \n"
      out$par <- paste0(
        "  // parameters for multivariate linear models \n",
        "  vector<lower=0>[K_trait] sigma; \n",
        "  cholesky_factor_corr[K_trait] Lrescor; \n")
      out$loop <- c(paste0(
        "  // restructure linear predictor \n",
        "  for (m in 1:N_trait) { \n",  
        "    for (k in 1:K_trait) { \n", 
        "      int n; \n",
        "      n <- (k - 1) * N_trait + m; \n"), 
        "    } \n  } \n")
      out$prior <- paste0(
        stan_prior(class = "sigma", coef = response, prior = prior),
        stan_prior(class = "Lrescor", prior = prior))
      if (family$family == "gaussian") {
        out$transD <- "  cholesky_factor_cov[K_trait] LSigma; \n"
        out$transC <- paste0(
          "  // compute cholesky factor of residual covariance matrix \n",
          "  LSigma <- diag_pre_multiply(sigma, Lrescor); \n")
      } else if (family$family %in% c("student", "cauchy")) {
        out$transD <- "  cov_matrix[K_trait] Sigma; \n"
        out$transC <- paste0(
          "  // compute residual covariance matrix \n",
          "  Sigma <- multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor)); \n")
      }
      out$genD <- paste0(
        "  matrix[K_trait,K_trait] Rescor; \n",
        "  vector<lower=-1,upper=1>[NC_trait] rescor; \n")
      out$genC <- paste0(
        "  // take only relevant parts of residual correlation matrix \n",
        "  Rescor <- multiply_lower_tri_self_transpose(Lrescor); \n",
        collapse(ulapply(2:nresp, function(i) lapply(1:(i-1), function(j)
          paste0("  rescor[",(i-1)*(i-2)/2+j,"] <- Rescor[",j,", ",i,"]; \n")))))
    } else if (!is.forked(family) && !is.categorical(family)) {
      stop("invalid multivariate model", call. = FALSE)
    }
  }
  out
}

stan_ordinal <- function(family, prior = prior_frame(), 
                         cse = FALSE, threshold = "flexible") {
  # Ordinal effects in Stan
  #
  # Args:
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   cse: logical; are there category specific effects?
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  stopifnot(is(family, "family"))
  out <- list()
  if (is.ordinal(family)) {
    # define Stan code similar for all ordinal models
    out$data <- "  int ncat;  // number of categories \n"
    th <- function(k, fam = family) {
      # helper function generating stan code inside ilink(.)
      sign <- ifelse(fam %in% c("cumulative", "sratio")," - ", " + ")
      ptl <- ifelse(cse, paste0(sign, "etap[k]"), "") 
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
        "    temp_Intercept[k] <- temp_Intercept1 + (k - 1.0) * delta; \n",
        "  } \n")
      out$prior <- paste0(stan_prior(class = "temp_Intercept1", prior = prior), 
                          stan_prior(class = "delta", prior = prior))
    }
    
    # generate Stan code specific for each ordinal model
    if (!(family == "cumulative" && ilink == "inv_logit")) {
      cse_arg <- ifelse(!cse, "", "row_vector etap, ")
      out$fun <- paste0(
        "  /* ", family, " log-PDF for a single response \n",
        "   * Args: \n",
        "   *   y: response category \n",
        "   *   eta: linear predictor \n",
        "   *   etap: optional linear predictor for category specific effects \n",
        "   *   thres: ordinal thresholds \n",
        "   * Returns: \n", 
        "   *   a scalar to be added to the log posterior \n",
        "   */ \n",
        "   real ", family, "_log(int y, real eta, ", cse_arg, "vector thres) { \n",
        "     int ncat; \n",
        "     vector[num_elements(thres) + 1] p; \n",
        if (family != "cumulative") "     vector[num_elements(thres)] q; \n",
        "     ncat <- num_elements(thres) + 1; \n")
      
      # define actual function content
      if (family == "cumulative") {
        out$fun <- paste0(out$fun,
        "     p[1] <- ", ilink, "(", th(1), "); \n",
        "     for (k in 2:(ncat - 1)) { \n", 
        "       p[k] <- ", ilink, "(", th("k"), ") - ",
        ilink, "(", th("k - 1"), "); \n", 
        "     } \n",
        "     p[ncat] <- 1 - ",ilink, "(", th("ncat - 1"), "); \n")
      } else if (family %in% c("sratio", "cratio")) {
        sc <- ifelse(family == "sratio", "1 - ", "")
        out$fun <- paste0(out$fun,
        "     for (k in 1:(ncat - 1)) { \n",
        "       q[k] <- ", sc, ilink, "(", th("k"), "); \n",
        "       p[k] <- 1 - q[k]; \n",
        "       for (kk in 1:(k - 1)) p[k] <- p[k] * q[kk]; \n", 
        "     } \n",
        "     p[ncat] <- prod(q); \n")
      } else if (family == "acat") {
        if (ilink == "inv_logit") {
          out$fun <- paste0(out$fun,
          "     p[1] <- 1.0; \n",
          "     for (k in 1:(ncat - 1)) { \n",
          "       q[k] <- ", th("k"), "; \n",
          "       p[k + 1] <- q[1]; \n",
          "       for (kk in 2:k) p[k + 1] <- p[k + 1] + q[kk]; \n",
          "       p[k + 1] <- exp(p[k + 1]); \n",
          "     } \n",
          "     p <- p / sum(p); \n")
        } else {
          out$fun <- paste0(out$fun,    
          "     for (k in 1:(ncat - 1)) \n",
          "       q[k] <- ", ilink, "(", th("k"), "); \n",
          "     for (k in 1:ncat) { \n",     
          "       p[k] <- 1.0; \n",
          "       for (kk in 1:(k - 1)) p[k] <- p[k] * q[kk]; \n",
          "       for (kk in k:(ncat - 1)) p[k] <- p[k] * (1 - q[kk]); \n",      
          "     } \n",
          "     p <- p / sum(p); \n")
        }
      }
      out$fun <- paste(out$fun, "    return categorical_log(y, p); \n   } \n")
    }
  }
  out
}

stan_categorical <- function(family) {
  # Stan code specific to categorical models
  # Args:
  #  family: the model family
  # Returns:
  #   a list of character strings defining the stan code
  #   specific for categorical models
  stopifnot(is(family, "family"))
  out <- list()
  if (is.categorical(family)) {
    out$data <- "  #include 'data_categorical.stan' \n" 
    out$tdataD <- "  vector[1] zero; \n"
    out$tdataC <- "  zero[1] <- 0; \n"
  }
  out
}

stan_zero_inflated_hurdle <- function(family) {
  # Stan code for zero-inflated and hurdle models
  # Args:
  #   family: the model family
  # Returns:
  #   a list of character strings defining the stan code
  #   specific for zero-inflated and hurdle models
  stopifnot(is(family, "family"))
  out <- list()
  if (is.zero_inflated(family) || is.hurdle(family)) {
    if (family$family == "zero_inflated_poisson") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_zero_inflated_poisson.stan' \n")
    } else if (family$family == "zero_inflated_negbinomial") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_zero_inflated_negbinomial.stan' \n")
    } else if (family$family == "zero_inflated_binomial") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_zero_inflated_binomial.stan' \n")
    } else if (family$family == "zero_inflated_beta") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_zero_inflated_beta.stan' \n")
    } else if (family$family == "hurdle_poisson") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_hurdle_poisson.stan' \n")
    } else if (family$family == "hurdle_negbinomial") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_hurdle_negbinomial.stan' \n")
    } else if (family$family == "hurdle_gamma") {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_hurdle_gamma.stan' \n")
    } 
  }
  out
}

stan_2PL <- function(family) {
  stopifnot(is(family, "family"))
  out <- list()
  if (is.2PL(family)) {
    out$transD <- "  vector[N_trait] eta_2PL;  // 2PL linear predictor \n"
    out$transC <- paste0("  eta_2PL <- head(eta, N_trait)", 
                         " .* exp(tail(eta, N_trait)); \n")
  }
  out
}

stan_inv_gaussian <- function(family, weights = FALSE, cens = FALSE, 
                              trunc = FALSE) {
  # stan code for inverse gaussian models
  #
  # Args:
  #   family: the model family
  #   weights: weights present?
  #   cens: censored data?
  #   trunc: truncated data?
  #
  # Returns:
  #   a list of character strings defining the stan code
  #   specific for inverse gaussian models
  stopifnot(is(family, "family"))
  out <- list()
  if (family$family == "inverse.gaussian") {
    out$data <- paste0(
      "  // quantities for the inverse gaussian distribution \n",
      "  vector[N] sqrt_Y;  // sqrt(Y) \n")
    if (weights || cens || trunc) {
      out$data <- paste0(out$data, "  vector[N] log_Y;  // log(Y) \n")
      out$fun <- paste0(out$fun, 
        "  #include 'fun_inv_gaussian_pw.stan' \n")
    } else {
      out$data <- paste0(out$data, "  real log_Y;  // sum(log(Y)) \n")
      out$fun <- paste0(out$fun, 
        "  #include 'fun_inv_gaussian_vector.stan' \n")
    } 
    if (cens || trunc) {
      out$fun <- paste0(out$fun, 
        "  #include 'fun_inv_gaussian_cdf.stan' \n",
        "  #include 'fun_inv_gaussian_ccdf.stan' \n")
    }
  }
  out
}

stan_disp <- function(disp, family = gaussian()) {
  # stan code for models with addition argument 'disp'
  # Args:
  #   disp: logical; are dispersion factors specified?
  #   family: the model family
  stopifnot(is(family, "family"))
  out <- list()
  if (disp) {
    par <- if (has_sigma(family)) "sigma"
           else if (has_shape(family)) "shape"
           else stop("invalid family for addition argument 'disp'")
    out$data <- "  vector<lower=0>[N] disp;  // dispersion factors \n"
    out$transD <- paste0("  vector<lower=0>[N] disp_", par, "; \n")
    out$transC <- paste0("  disp_", par, " <- ", par, " * disp; \n")
  }
  out
}

stan_misc_functions <- function(family = gaussian(), kronecker = FALSE) {
  # stan code for user defined functions
  #
  # Args:
  #   family: the model family
  #   kronecker: logical; is the kronecker product needed?
  #
  # Returns:
  #   a string containing defined functions in stan code
  stopifnot(is(family, "family"))
  out <- NULL
  if (family$link == "cauchit") {
    out <- paste0(out, "  #include 'fun_cauchit.stan' \n")
  }
  if (kronecker) {
    out <- paste0(out, "  #include 'fun_as_matrix.stan' \n",
                  "  #include 'fun_kronecker.stan' \n")
  }
  out
}

stan_prior <- function(class, coef = NULL, group = NULL, nlpar = NULL, 
                       suffix = "", matrix = FALSE, prior = prior_frame(), 
                       wsp = 2) {
  # Define priors for parameters in Stan language
  # 
  # Args:
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   gi: index of the grouping factor
  #   nlpar: the name of a non-linear parameter
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   matrix: logical; corresponds the class to a parameter matrix?
  #   wsp: an integer >= 0 defining the number of spaces 
  #      in front of the output string
  # 
  # Returns:
  #   A character strings in stan language that defines priors for a given class of parameters
  #   If a parameter has has no corresponding prior in prior 
  #   and also no internal default in stan_prior, an empty string is returned.
  
  # only consider user defined priors related to this class and group
  wsp <- collapse(rep(" ", wsp))
  prior_only <- isTRUE(attr(prior, "prior_only"))
  keep <- which(prior$class == class & (prior$coef %in% coef | !nchar(prior$coef)))
  user_prior <- prior[keep, ]
  if (!is.null(group)) {
    keep_group <- which(user_prior$group == group | !nchar(user_prior$group))
    user_prior <- user_prior[keep_group, ]
  }
  if (!is.null(nlpar)) {
    keep_nlpar <- which(user_prior$nlpar == nlpar | !nchar(user_prior$nlpar))
    user_prior <- user_prior[keep_nlpar, ]
  }
  if (!nchar(class) && nrow(user_prior)) {
    # increment_log_prob statements are directly put into the Stan code
    return(collapse(wsp, user_prior$prior, "; \n"))
  } 
  
  # get base prior
  igroup <- which(with(user_prior, !nchar(coef) & nchar(group) & nchar(prior)))
  inlpar <- which(with(user_prior, !nchar(coef) & nchar(nlpar) & nchar(prior)))
  iclass <- which(with(user_prior, !nchar(coef) & !nchar(group) & nchar(prior)))
  if (length(igroup)) {  
    # if there is a global prior for this group
    base_prior <- user_prior[igroup, "prior"]
  } else if (length(inlpar)) {
    base_prior <- user_prior[inlpar, "prior"]
  } else if (length(iclass)) {  
    # if there is a global prior for this class
    base_prior <- user_prior[iclass, "prior"]
  } else {  
    # no proper prior for this class
    base_prior <- ""
  } 
  
  individual_prior <- function(i, max_index) {
    # individual priors for each parameter of a class
    if (max_index > 1L || matrix) {
      index <- paste0("[",i,"]")      
    } else {
      index <- ""
    }
    uc_prior <- user_prior$prior[match(coef[i], user_prior$coef)]
    if (!is.na(uc_prior) & nchar(uc_prior)) { 
      # user defined prior for this parameter
      coef_prior <- uc_prior
    } else {
      # base prior for this parameter
      coef_prior <- base_prior  
    }  
    if (nchar(coef_prior) > 0) {  # implies a proper prior
      out <- paste0(wsp, class, index, " ~ ", coef_prior, "; \n")
    } else {
      out <- "" # implies an improper flat prior
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(class, suffix)
  if (any(with(user_prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(1:length(coef), individual_prior, max_index = length(coef))
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- paste0("to_vector(", class, ")")
    }
    out <- paste0(wsp, class, " ~ ", base_prior, "; \n")
  } else {
    out <- ""
  }
  if (class == "b" && !is.null(attr(prior, "hs_df"))) {
    # add horseshoe shrinkage priors
    hs_shrinkage_priors <- paste0(
      "  hs_local ~ student_t(", attr(prior, "hs_df"), ", 0, 1); \n",
      "  hs_global ~ cauchy(0, 1); \n")
    out <- c(hs_shrinkage_priors, out)
  }
  out <- collapse(out)
  if (prior_only && nchar(class) && !nchar(out)) {
    stop(paste0("Sampling from priors is not possible because ", 
                "not all parameters have proper priors. \n",
                "Error occured for class '", class, "'."), 
         call. = FALSE)
  }
  out
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
    take <- !grepl("^(z|temp)_|^increment_log_prob\\(", pars)
    pars <- rename(pars[take], symbols = c("^L_", "^Lrescor"), 
                   subs = c("cor_", "rescor"), fixed = FALSE)
    dis <- gsub("~", "", regmatches(prior, gregexpr("~[^\\(]+", prior))[[1]])[take]
    args <- regmatches(prior, gregexpr("\\([^;~]+\\);", prior))[[1]][take]
    type <- rep("real", length(pars))
    
    # rename parameters containing indices
    has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
    pars[has_ind] <- sapply(pars[has_ind], function(par) {
      ind <- regmatches(par, gregexpr("\\[[[:digit:]]+\\]", par))
      ind <- as.numeric(substr(ind, 2, nchar(ind) - 1))
      if (grepl("^b\\[", par)) {
        par <- paste0("b_",ind)
      } else if (grepl("^bp\\[", par)) {
        par <- paste0("bp_",ind)
      } else if (grepl("^sigma\\[", par)) {
        par <- paste0("sigma_",ind)
      } else if (grepl("^sd_", par)) {
        par <- gsub("\\[[[:digit:]]+\\]", paste0("_",ind), par)
      }
      return(par)
    })
    
    # special treatment of lkj_corr_cholesky priors
    args <- ifelse(grepl("corr_cholesky$", dis), 
                   paste0("(2,", substr(args, 2, nchar(args)-1), "[1,2];"), 
                   args)
    dis <- sub("corr_cholesky$", "corr", dis)
    
    # special treatment of simplex parameters
    which_simplex <- which(grepl("^simplex_", pars))
    for (i in seq_along(which_simplex)) {
      type[which_simplex[i]] <- paste0("vector[rows(con_simplex_", i, ")]")
    }
    
    # extract possible boundaries
    par_declars <- unlist(strsplit(par_declars, "\n", fixed = TRUE))
    all_pars <- get_matches(" [^[:blank:]]+;", par_declars) 
    all_pars <- substr(all_pars, 2, nchar(all_pars) - 1)
    all_bounds <- get_matches("<.+>", par_declars, simplify = FALSE)
    all_bounds <- ulapply(all_bounds, function(x) if (length(x)) x else "")
    bounds <- rep("", length(pars))
    for (i in seq_along(all_pars)) {
      k <- which(grepl(paste0("^", all_pars[i]), pars))
      bounds[k] <- all_bounds[i]
    }
    has_bounds <- as.logical(nchar(bounds))
    
    # distinguish between bounded and unbounded parameters
    #bound <- grepl("^sd|^sigma|^shape$|^nu$|^hs_local$|^hs_global$", pars) |  
    #  family$family == "cumulative" & grepl("^delta$", pars)
    if (any(has_bounds)) {  
      # bounded parameters have to be sampled in the model block
      #lower_bound <- ifelse(pars[bound] == "nu", 1, 0)
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
        args[match("b", pars)] <- "(0, prior_hs_local * prior_hs_global);" 
      } 
      out$genD <- collapse(
        "  ", type[no_bounds], " prior_", pars[no_bounds], "; \n")
      out$genC <- paste0(
        "  // additionally draw samples from priors \n",
        collapse("  prior_", pars[no_bounds], " <- ",
                 dis[no_bounds], "_rng", args[no_bounds], " \n"))
    }
  }
  out
}

stan_ilink <- function(link) {
  # find the inverse link to a given link function
  # 
  # Args:
  #   link: the link function
  #
  # Returns: 
  #   the inverse link function for stan; a character string
  switch(link, identity = "", log = "exp", inverse = "inv", 
         sqrt = "square", "1/mu^2" = "inv_sqrt",
         logit = "inv_logit", probit = "Phi", 
         probit_approx = "Phi_approx", cloglog = "inv_cloglog", 
         cauchit = "inv_cauchit")
}

stan_eta_fixef <- function(fixef, sparse = FALSE) {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the fixed effects
  #   sparse: logical; use sparse matrix multiplication?
  if (length(fixef)) {
    if (sparse) {
      eta_fixef <- "csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b)"
    } else {
      eta_fixef <- "X * b"
    }
  } else { 
    eta_fixef <- "rep_vector(0, N)"
  }
  eta_fixef
}

stan_eta_ranef <- function(ranef, nlpar = NULL) {
  # Write the random effects part of the linear predictor
  # Args:
  #   ranef: a named list returned by gather_ranef
  #   par: an optional character string to add to the variable names
  #        (used for non-linear models)
  # Returns:
  #   A string containing the random effects part of the linear predictor
  eta_ranef <- ""
  for (i in seq_along(ranef)) {
    nli <- if (!is.null(nlpar)) paste0(nlpar, "_", i) else i
    if (length(ranef[[i]]) == 1L) {
      eta_ranef <- paste0(eta_ranef, 
        " + r_", nli,"[J_", nli,"[n]] * Z_", nli,"[n]")
    } else {
      j <- seq_along(ranef[[i]])
      eta_ranef <- paste0(eta_ranef, collapse(
        " + r_", nli, "_", j, "[J_", nli,"[n]]",
        " * Z_", nli, "_", j, "[n]"))
    }
  }
  eta_ranef
}

stan_eta_monef <- function(monef) {
  eta_monef <- ""
  for (i in seq_along(monef)) {
    eta_monef <- paste0(eta_monef,
      " + bm[", i, "] * monotonous(simplex_", i, ", Xm[n, ", i, "])")
  }
  eta_monef
}

stan_eta_transform <- function(family, link, add = FALSE) {
  # indicate whether eta needs to be transformed
  # in the transformed parameters block
  # Args:
  #   add: is the model weighted, censored, or truncated?
  !(add || !is.skewed(family) && link == "identity" 
    || family %in% "gaussian" && link == "log"
    || is.count(family) && link == "log" 
    || is.binary(family) && link == "logit"
    || is.ordinal(family) || is.categorical(family) 
    || is.zero_inflated(family) || is.hurdle(family))
}

stan_eta_ilink <- function(family, link, disp = FALSE) {
  # correctly apply inverse link to eta
  ilink <- stan_ilink(link)
  shape <- ifelse(disp, "disp_shape[n]", "shape")
  fl <- ifelse(family %in% c("gamma", "exponential"), 
               paste0(family,"_",link), family)
  switch(fl, c(paste0(ilink,"("), ")"),
         gamma_log = c(paste0(shape, " * exp(-("), "))"),
         gamma_inverse = c(paste0(shape, " * ("), ")"),
         gamma_identity = c(paste0(shape, " / ("), ")"),
         exponential_log = c("exp(-(", "))"),
         exponential_inverse = c("(", ")"),
         exponential_identity = c("inv(", ")"),
         weibull = c(paste0(ilink,"(("), 
                     paste0(") / ", shape, ")")))
}

has_build_in_fun <- function(family, link) {
  # indicates if a family-link combination has a build in 
  # function in Stan (such as binomial_logit)
  (family %in% c("binomial", "bernoulli", "cumulative", "categorical")
   && link == "logit" || is.count(family) && link == "log")
}

needs_kronecker <- function(ranef, names_cov_ranef) {
  # checks if a model needs the kronecker product
  # Args: 
  #   ranef: named list returned by gather_ranef
  #   names_cov_ranef: names of the grouping factors that
  #                    have a cov.ranef matrix 
  .fun <- function(x, names) {
    length(x) > 1 && attr(x, "group") %in% names && attr(x, "cor")
  }
  any(sapply(ranef, .fun, names = names_cov_ranef))
}

