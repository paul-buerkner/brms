stan_fixef <- function(fixef, paref, family = gaussian(), 
                       prior = prior_frame(), has_intercept = TRUE, 
                       threshold = "flexible") {
  # Stan code for fixec effects
  #
  # Args:
  #   fixef: names of the fixed effects
  #   paref: names of the category specific effects
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   has_intercept: logical; fixed effects intercept present?
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   a list containing Stan code related to fixed effects
  out <- list()
  if (has_intercept) {
    if (is.categorical(family)) {
      out$par <- paste0(out$par,
                        "  row_vector [ncat - 1] temp_Intercept;",
                        "  // temporary intercepts \n")
      out$genD <- paste0("  row_vector[ncat - 1] b_Intercept;",
                         "  // fixed effects intercepts \n")
      subtract <- ifelse(length(paref), " - to_row_vector(Xp_means) * bp", "")
      out$genC <- paste0("  b_Intercept <- temp_Intercept", subtract, "; \n")
      out$prior <- stan_prior("temp_Intercept", prior = prior)
    } else if (is.ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      subtract <- ifelse(length(fixef), " - dot_product(X_means, b)", "") 
      out$genC <- paste0("  b_Intercept <- temp_Intercept", subtract, "; \n")
    } else {
      out$par <- paste0(out$par,
                        "  real temp_Intercept;  // temporary Intercept \n")
      out$genD <- "  real b_Intercept;  // fixed effects intercept \n"
      subtract <- ifelse(length(fixef), " - dot_product(X_means, b)", "") 
      out$genC <- paste0("  b_Intercept <- temp_Intercept", subtract, "; \n")
      out$prior <- stan_prior("temp_Intercept", prior = prior)
    }
  }
  if (length(fixef)) {
    out$data <- paste0(out$data, 
      "  int<lower=1> K;  // number of fixed effects \n", 
      "  matrix[N, K] X;  // FE design matrix \n",
      "  vector[K] X_means;  // column means of X \n")
    bound <- with(prior, bound[class == "b" & coef == ""])
    out$par <- paste0(out$par,
      "  vector", bound, "[K] b;  // fixed effects \n") 
    fixef_prior <- stan_prior(class = "b", coef = fixef, prior = prior)
    out$prior <- paste0(out$prior, fixef_prior)
  }
  if (length(paref)) {
    out$data <- paste0(out$data, 
      "  int<lower=1> Kp;  // number of category specific effects \n",
      "  matrix[N, Kp] Xp;  // CSE design matrix \n",
      if (is.categorical(family)) 
      "  vector[Kp] Xp_means;  // column means of Xp \n")
    bound <- with(prior, bound[class == "b" & coef == ""])
    out$par <- paste0(out$par,
      "  matrix", bound, "[Kp, ncat - 1] bp;  // category specific effects \n")
    paref_prior <- stan_prior(class = "bp", coef = paref, prior = prior)
    out$prior <- paste0(out$prior, paref_prior)
  }
  out
}

stan_ranef <- function(i, ranef, prior = prior_frame(), 
                       names_cov_ranef = NULL, par = "") {
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
  pi <- if (nchar(par)) paste0(par, "_", i) else i
  out <- list()
  out$data <- paste0(
    "  // data for random effects of ", g, " \n",
    "  int<lower=1> J_", pi, "[N];  // RE levels \n",
    "  int<lower=1> N_", pi, ";  // number of levels \n",
    "  int<lower=1> K_", pi, ";  // number of REs \n",
    if (ccov) paste0(
      "  matrix[N_", pi, ", N_", pi,"] cov_", pi,";",
      "  // user defined covariance matrix \n"))
  
  out$prior <- stan_prior(class = "sd", group = pi, coef = r, prior = prior)
  if (length(r) == 1) {  # only one random effect
    out$data <- paste0(out$data, "  real Z_", pi, "[N];  // RE design matrix \n")
    out$par <- paste0("  vector[N_", pi, "] pre_", pi, ";  // unscaled REs \n",
                      "  real<lower=0> sd_", pi, ";  // RE standard deviation \n")
    out$prior <- paste0(out$prior,"  pre_", pi, " ~ normal(0, 1); \n")
    out$transD <- paste0("  vector[N_", pi, "] r_", pi, ";  // REs \n")
    out$transC <- paste0("  r_", pi,  " <- sd_", pi, " * (", 
                         if (ccov) paste0("cov_", pi, " * "), "pre_", pi, ");",
                         "  // scale REs \n")
  } else if (length(r) > 1 && cor) {  
    # multiple correlated random effects
    out$data <- paste0(out$data,  
      "  row_vector[K_", pi, "] Z_", pi, "[N];  // RE design matrix \n",  
      "  int NC_", pi, ";  // number of correlations \n")
    out$par <- paste0(
      "  matrix[N_", pi, ", K_", pi, "] pre_", pi, ";  // unscaled REs \n",
      "  vector<lower=0>[K_", pi, "] sd_", pi, ";  // RE standard deviation \n",
      "  // cholesky factor of correlation matrix \n",
      "  cholesky_factor_corr[K_", pi, "] L_", pi, ";")
    out$prior <- paste0(out$prior, 
                        stan_prior(class = "L", group = pi,  prior = prior),
                        "  to_vector(pre_", pi, ") ~ normal(0, 1); \n")
    out$transD <- paste0("  vector[K_", pi, "] r_", pi, "[N_", pi, "];  // REs \n")
    if (ccov) {  # customized covariance matrix supplied
      out$transC <- paste0("  r_", pi, 
        " <- to_array(kronecker_cholesky(cov_", pi, ", L_", pi, ", sd_", pi, ") * ",
        "to_vector(pre_", pi, "), N_", pi, ", K_", pi, ");  // scale REs \n")
    } else { 
      out$transC <- paste0(
        "  for (i in 1:N_", pi, ") { \n",
        "    r_", pi,  "[i] <- sd_", pi, " .* (L_", pi, " * ", 
        "to_vector(pre_", pi, "[i]));  // scale REs \n  } \n")
    }
    # return correlations above the diagonal only
    cors_genC <- ulapply(2:length(r), function(k) 
      lapply(1:(k - 1), function(j) paste0(
        "  cor_", pi, "[", (k - 1) * (k - 2) / 2 + j, 
        "] <- Cor_", pi, "[", j, ",", k, "]; \n")))
    out$genD <- paste0(
      "  corr_matrix[K_", pi, "] Cor_", pi, "; \n",
      "  vector<lower=-1,upper=1>[NC_", pi, "] cor_", pi, "; \n")
    out$genC <- paste0(
      "  // take only relevant parts of correlation matrix \n",
      "  Cor_", pi, " <- multiply_lower_tri_self_transpose(L_", pi, "); \n",
      collapse(cors_genC)) 
  } else if (length(r) > 1 && !cor) {
    # multiple uncorrelated random effects
    j <- seq_along(r)
    out$data <- paste0(out$data, "  matrix[N, K_", pi, "] Z_", pi, ";",
                       "  // RE design matrix \n")
    out$par <- paste0("  vector[N_", pi, "] pre_", pi, "[K_", pi, "];",
                      "  // unscaled REs \n",
                      "  vector<lower=0>[K_", pi, "] sd_", pi, ";",
                      "  // RE standard deviation \n")
    out$prior <- paste0(out$prior, 
                        collapse("  pre_", pi, "[", j, "] ~ normal(0, 1); \n"))
    out$transD <- collapse("  vector[N_", pi, "] r_", pi, "_", j, ";  // REs \n")
    out$transC <- collapse(
      "  r_", pi, "_",j," <- sd_", pi, "[", j, "] * (", 
      if (ccov) paste0("cov_", pi, " * "), "pre_", pi, "[", j, "]);",
      "  // scale REs \n")
  }
  out
}

stan_llh <- function(family, se = FALSE, weights = FALSE, trials = FALSE, 
                     cens = FALSE, disp = FALSE, trunc = .trunc(), 
                     autocor = cor_arma(), partial = FALSE, is_multi = FALSE) {
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
  #   partial: a flag whether catgory specific effects are present
  #            (for ordinal models only)
  #   is_multi: is the model multivariate?
  #
  # Returns:
  #   a string containing the likelihood of the model in stan language
  if (!is(family, "family"))
    stop("family must be of class family")
  link <- family$link
  type <- family$type
  family <- family$family
  is_linear <- is.linear(family)
  is_catordinal <- is.ordinal(family) || is.categorical(family)
  is_count <- is.count(family)
  is_skewed <- is.skewed(family)
  is_binary <- is.binary(family)
  is_hurdle <- is.hurdle(family)
  is_zero_inflated <- is.zero_inflated(family)
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
  if (!is.null(type)) family <- paste0(family,"_",type)
  
  simplify <- !is_trunc && !cens && 
    (is_binary && link == "logit" || is_count && link == "log" ||
       family %in% c("cumulative", "categorical") && link == "logit") 
  n <- ifelse(cens || weights || is_trunc || is_catordinal ||
                is_hurdle || is_zero_inflated, "[n]", "")
  ns <- ifelse((se || trials || disp) && (cens || weights || is_trunc) 
               || (trials && is_zero_inflated), "[n]", "")
  disp <- ifelse(disp, "disp_", "")
  sigma <- paste0(ifelse(se, "se", paste0(disp, "sigma")), ns)
  shape <- paste0(disp, "shape", ns)
  ordinal_args <- paste0("eta[n], ", if (partial) "etap[n], ", 
                         "temp_Intercept")
  # use inverse link in likelihood statement only 
  # if it does not prevent vectorization 
  ilink <- ifelse(n == "[n]" && !simplify, stan_ilink(link), "")
  if (n == "[n]") {
    if (is_hurdle || is_zero_inflated) {
      eta <- paste0(ilink,"(eta[n]), ",ilink,"(eta[n + N_trait])")
    } else {
      fl <- ifelse(family %in% c("gamma", "exponential"), 
                   paste0(family,"_",link), family)
      eta <- switch(fl, paste0(ilink,"(eta[n])"),
                    gamma_log = paste(shape, "* exp(-eta[n])"),
                    gamma_inverse = paste(shape, "* eta[n]"),
                    gamma_identity =  paste(shape, "/ eta[n]"),
                    exponential_log = "exp(-eta[n])",
                    exponential_inverse = "eta[n]",
                    exponential_identity = "inv(eta[n])",
                    weibull = paste0(ilink, "(eta[n] / ", shape, ")"),
                    bernoulli_2PL = paste0(ilink, "(eta_2PL[n])"))
    }
  } else {
    # possible transformations already performed
    # in the transformed parameters block
    eta <- "eta"
  }
  
  if (simplify) { 
    llh_pre <- switch(family,
      poisson = c("poisson_log", paste0("eta",n)), 
      negbinomial = c("neg_binomial_2_log", paste0("eta",n,", shape")),
      geometric = c("neg_binomial_2_log", paste0("eta",n,", 1")),
      cumulative = c("ordered_logistic", "eta[n], temp_Intercept"),
      categorical = c("categorical_logit", 
                      "to_vector(append_col(zero, eta[n] + etap[n]))"), 
      binomial = c("binomial_logit", paste0("trials",ns,", eta",n)), 
      bernoulli = c("bernoulli_logit", paste0("eta",n)),
      bernoulli_2PL = c("bernoulli_logit", paste0("eta_2PL",n)))
  } else {
    llh_pre <- switch(family,
      gaussian = c("normal", paste0(eta, ", ", sigma)),
      gaussian_cov = c("normal_cov", paste0(eta,", se2, N_tg, ", 
                       "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      student = c("student_t",  paste0("nu, ",eta, ", ", sigma)),
      student_cov = c("student_t_cov", paste0("nu, ",eta,", se2, N_tg, ", 
                      "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      cauchy = c("cauchy", paste0(eta,", ", sigma)),
      cauchy_cov = c("student_t_cov", paste0("1, ",eta,", se2, N_tg, ", 
                     "begin_tg, end_tg, nobs_tg, res_cov_matrix")),
      lognormal = c("lognormal", paste0(eta,", sigma",ns)),
      multi_gaussian = c("multi_normal_cholesky", paste0("Eta",n,", LSigma")),
      multi_student = c("multi_student_t", paste0("nu, Eta",n,", Sigma")),
      multi_cauchy = c("multi_student_t", paste0("1.0, Eta",n,", Sigma")),
      poisson = c("poisson", eta),
      negbinomial = c("neg_binomial_2", paste0(eta, ", ", shape)),
      geometric = c("neg_binomial_2", paste0(eta,", 1")),
      binomial = c("binomial", paste0("trials",ns,", ",eta)),
      bernoulli = c("bernoulli", eta), 
      bernoulli_2PL = c("bernoulli", eta), 
      gamma = c("gamma", paste0(shape, ", ", eta)), 
      exponential = c("exponential", eta),
      weibull = c("weibull", paste0(shape, ", ", eta)), 
      inverse.gaussian = c(paste0("inv_gaussian", if (!nchar(n)) "_vector"), 
                           paste0(eta, ", shape, log_Y",n,", sqrt_Y",n)),
      beta = c("beta", paste0(eta, " * phi, (1 - ", eta, ") * phi")),
      categorical = c("categorical", "p[n]"),
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
      ifelse(!weights, paste0("Y[n] ~ ", llh_pre[1],"(",llh_pre[2],"); \n"),
             paste0("increment_log_prob(", add_weights, 
                    llh_pre[1], "_log(Y[n], ",llh_pre[2],")); \n")),
      "    else { \n",         
      "      if (cens[n] == 1) increment_log_prob(", add_weights, 
      llh_pre[1], "_ccdf_log(Y[n], ",llh_pre[2],")); \n",
      "      else increment_log_prob(", add_weights, 
      llh_pre[1], "_cdf_log(Y[n], ",llh_pre[2],")); \n",
      "    } \n"),
    weights = paste0("  lp_pre[n] <- ", llh_pre[1], "_log(Y[n], ",
                     llh_pre[2],"); \n"),
    general = paste0("  Y", n, " ~ ", llh_pre[1],"(",llh_pre[2],")", 
                     code_trunc, "; \n")) 
  llh
}

stan_eta <- function(family, fixef, ranef = list(), paref = NULL, 
                     has_intercept = TRUE, autocor = cor_arma(),  
                     add = FALSE, offset = FALSE, is_multi = FALSE) {
  # linear predictor in Stan
  #
  # Args:
  #   family: the model family
  #   fixef: names of the fixed effects
  #   ranef: a named list returned by gather_ranef
  #   paref: names of the category specific effects
  #   has_intercept: has the model a fixed effects intercept?
  #   autocor: autocorrelation structure
  #   add: is the model weighted, censored, or truncated?
  #   offset: is an offset defined?
  #   is_multi: is the model multivariate?
  # 
  # Return:
  #   the linear predictor in stan language
  if (!is(family, "family"))
    stop("family must be of class family")
  link <- family$link
  family <- family$family
  is_ordinal <- is.ordinal(family)
  is_cat <- is.categorical(family)
  # initialize eta
  eta <- list()
  eta$transD <- paste0(
    "  vector[N] eta;  // linear predictor \n", 
    if (length(paref) || is.categorical(family)) 
      paste0("  matrix[N, ncat - 1] etap;",
             "  // linear predictor for category specific effects \n"),
    if (is_multi) 
      paste0("  vector[K_trait] Eta[N_trait];",
             "  // multivariate linear predictor matrix \n"))
  eta_obj <- ifelse(is_multi, "Eta[m, k]", "eta[n]")
  s <- ifelse(is_multi, "  ", "")
  
  # transform eta before it is passed to the likelihood
  eta$transform <- stan_eta_transform(family, link, add = add)
  eta_ilink <- rep("", 2)
  if (eta$transform || (get_ar(autocor) && !use_cov(autocor))) {
    eta_ilink <- stan_eta_ilink(family, link)
    if (get_ar(autocor)) {
      eta_ar <- ifelse(!use_cov(autocor), " + head(E[n], Kar) * ar", "")
      eta$transC3 <- paste0("    ", s, eta_obj," <- ", eta_ilink[1], 
                            eta_obj, eta_ar, eta_ilink[2], "; \n")
      # don't apply link function twice
      eta_ilink <- rep("", 2)
    }
  }
  
  # define fixed, random, and autocorrelation effects
  eta_re <- stan_eta_re(ranef)
  etap <- if (length(paref) || is_cat) {
    paste0("  etap <- ", 
           ifelse(length(paref), "Xp * bp", "rep_matrix(0, N, ncat - 1)"),
           if (is_cat && has_intercept) " + rep_matrix(temp_Intercept, N)", "; \n")
  }
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor), 
                   " + head(E[n], Kma) * ma", "")
  if (nchar(eta_re) || nchar(eta_ma) || is_multi || nchar(eta_ilink[1])) {
    eta$transC2 <- paste0("    ",s, eta_obj," <- ", eta_ilink[1], 
                          "eta[n]", eta_ma, eta_re, eta_ilink[2],"; \n")
    eta_ilink <- rep("", 2)
  }
  eta$transC1 <- paste0(
    "  // compute linear predictor \n",
    "  eta <- ", ifelse(length(fixef), "X * b", "rep_vector(0, N)"), 
    if (has_intercept && !(is_ordinal || is_cat)) " + temp_Intercept",
    if (offset) " + offset",
    if (get_arr(autocor)) " + Yarr * arr", 
    "; \n", etap)
  eta
}

stan_nonlinear <- function(effects, data, family = gaussian(), 
                           add = FALSE, cov_ranef = NULL, 
                           prior = prior_frame()) {
  # prepare Stan code for non-linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   add: Is the model weighted, censored, or truncated?
  #   cov_ranef: a list of user-defined covariance matrices
  #   prior: a prior_frame object
  out <- list()
  if (length(effects$nonlinear)) {
    out$data <- "  // data for non-linear fixed effects \n"
    out$par <- "  // non-linear fixed effects \n"
    for (i in seq_along(effects$nonlinear)) {
      nlp <- names(effects$nonlinear)[i]
      eta <- paste0("eta_", nlp)
      out$transD <- paste0(out$transD, "  vector[N] ", eta, "; \n")
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
          stan_prior(class = "b", coef = fixef, nlpar = nlp, prior = prior))
      } else {
        out$transC1 <- paste0(out$transC1, 
          "  ", eta, " <- rep_vector(0, N); \n")  
      }
      # include random effects
      ranef <- gather_ranef(effects$nonlinear[[i]], data = data)
      if (length(ranef)) {
        text_ranef <- lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
                             names_cov_ranef = names(cov_ranef), 
                             prior = prior, par = nlp)
        text_ranef <- collapse_lists(text_ranef)
        out$data <- paste0(out$data, text_ranef$data)
        out$prior <- paste0(out$prior, text_ranef$prior)
        out$par <- paste0(out$par, text_ranef$par)
        out$transD <- paste0(out$transD, text_ranef$transD)
        out$transC1 <- paste0(out$transC1, text_ranef$transC)
        out$transC2 <- paste0(out$transC2, 
          "    ", eta, "[n] <- ", eta, "[n]", 
          stan_eta_re(ranef, par = nlp), "; \n") 
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
    nlmodel <- gsub(" ", "", deparse(effects$fixed[[3]]))
    nlmodel <- wsp(rename(nlmodel, meta_sym, wsp(meta_sym))) 
    nlmodel <- rename(nlmodel, c(nlpars, covars, " ( ", " ) "), 
                      c(new_nlpars, new_covars, "(", ")"))
    # possibly transform eta in the transformed params block
    transform <- stan_eta_transform(family$family, family$link, add = add)
    if (transform) {
      eta_ilink <- stan_eta_ilink(family$family, family$link)
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
  if (!is(family, "family"))
    stop("family must be of class family")
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
      s <- ifelse(is_multi, "      ", "    ")
      link <- c(identity = "", log = "log", inverse = "inv")[family$link]
      out$transD <- paste0("  matrix[N, Karma] E;  // ARMA design matrix \n",
                           "  vector[N] e;  // residuals \n") 
      out$transC1 <- "  E <- E_pre; \n" 
      out$transC2 <- paste0(
        s, "// calculation of ARMA effects \n",
        s, "e[n] <- ", link, "(Y[", index, "]) - eta[n]", "; \n",
        s, "for (i in 1:Karma) { \n", 
        s, "  if (n + 1 - i > 0 && n < N && tg[n + 1] == tg[n + 1 - i]) { \n",
        s, "     E[n + 1, i] <- e[n + 1 - i]; \n",
        s, "  } \n",
        s, "} \n")
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
  if (!is(family, "family"))
    stop("family must be of class family")
  out <- list()
  nresp <- length(response)
  if (nresp > 1) {
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
    } else if (!is.forked(family)) {
      stop("invalid multivariate model", call. = FALSE)
    }
  }
  out
}

stan_ordinal <- function(family, prior = prior_frame(), 
                         partial = FALSE, threshold = "flexible") {
  # Ordinal effects in Stan
  #
  # Args:
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   partial: logical; are there partial effects?
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  if (!is(family, "family"))
    stop("family must be of class family")
  out <- list()
  if (is.ordinal(family)) {
    # define Stan code similar for all ordinal models
    th <- function(k, fam = family) {
      # helper function generating stan code inside ilink(.)
      sign <- ifelse(fam %in% c("cumulative", "sratio")," - ", " + ")
      ptl <- ifelse(partial, paste0(sign, "etap[k]"), "") 
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
      cse_arg <- ifelse(!partial, "", "row_vector etap, ")
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

stan_zero_inflated_hurdle <- function(family) {
  # stan code for zero-inflated and hurdle models
  #
  # Args:
  #   family: the model family
  #
  # Returns:
  #   a list of character strings defining the stan code
  #   specific for zero-inflated and hurdle models
  if (!is(family, "family"))
    stop("family must be of class family")
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
  if (!is(family, "family"))
    stop("family must be of class family")
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
  if (!is(family, "family"))
    stop("family must be of class family")
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
  out <- list()
  if (disp) {
    par <- if (has_sigma(family)) "sigma"
           else if (has_shape(family)) "shape"
           else stop("invalid family for addition argument 'disp'")
    out$data <- "  vector<lower=0>[N] disp;  // dispersion factors \n"
    out$transD <- paste0("  vector<lower=0>[N] disp_", par, ";")
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
  if (!is(family, "family"))
    stop("family must be of class family")
  out <- NULL
  if (family$link == "cauchit") {
    out <- paste0(out, "  #include 'fun_cauchit.stan' \n")
  }
  if (kronecker) {
    out <- paste0(out, "  #include 'fun_to_array.stan' \n",
                  "  #include 'fun_kronecker_cholesky.stan' \n")
  }
  out
}

stan_prior <- function(class, coef = NULL, group = NULL, nlpar = NULL,
                       prior = prior_frame(), s = 2) {
  # Define priors for parameters in Stan language
  # 
  # Args:
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   nlpar: the name of a non-linear parameter
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   s: an integer >= 0 defining the number of spaces 
  #      in front of the output string
  # 
  # Returns:
  #   A character strings in stan language that defines priors for a given class of parameters
  #   If a parameter has has no corresponding prior in prior 
  #   and also no internal default in stan_prior, an empty string is returned.
  
  # only consider user defined priors related to this class and group
  s <- collapse(rep(" ", s))
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
    return(collapse(s, user_prior$prior, "; \n"))
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
    if (max_index > 1 || class == "bp") {
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
      return(paste0(s, class, index, " ~ ", coef_prior, "; \n"))
    } else {
      return("")  # implies an improper flat prior
    }
  }
  
  if (!is.null(nlpar))
    class <- paste0(class, "_", nlpar)
  if (!is.null(group)) 
    class <- paste0(class, "_", group)
  # generate stan prior statements
  if (any(with(user_prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(1:length(coef), individual_prior, max_index = length(coef))
  } else if (nchar(base_prior) > 0) {
    if (class == "bp") {
      class <- "to_vector(bp)"
    }
    out <- paste0(s, class, " ~ ", base_prior, "; \n")
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
  return(collapse(out))
}

stan_rngprior <- function(sample_prior, prior, family = gaussian(),
                          hs_df = NULL) {
  # stan code to sample from priors seperately
  #
  # Args:
  #   sample_prior: take samples from priors?
  #   prior: the character string taken from stan_prior
  #   family: the model family
  #   hs_df: hs_df degrees of freedom
  #
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  if (!is(family, "family"))
    stop("family must be of class family")
  out <- list()
  if (sample_prior) {
    prior <- gsub(" ", "", paste0("\n",prior))
    pars <- gsub("\\\n|to_vector\\(|\\)", "", 
                 regmatches(prior, gregexpr("\\\n[^~]+", prior))[[1]])
    take <- !grepl("^pre_|^increment_log_prob\\(", pars)
    pars <- rename(pars[take], symbols = c("^L_", "^Lrescor"), 
                   subs = c("cor_", "rescor"), 
                   fixed = FALSE)
    dis <- gsub("~", "", regmatches(prior, gregexpr("~[^\\(]+", prior))[[1]])[take]
    args <- regmatches(prior, gregexpr("\\([^;~]+\\);", prior))[[1]][take]
    
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
    
    # distinguish between bounded and unbounded parameters
    # do not change | to ||
    bound <- grepl("^sd|^sigma|^shape$|^nu$|^hs_local$|^hs_global$", pars) |  
      family$family == "cumulative" & grepl("^delta$", pars)
    if (any(bound)) {  
      # bounded parameters have to be sampled in the model block
      lower_bound <- ifelse(pars[bound] == "nu", 1, 0)
      out$par <- paste0("  // parameters to store prior samples \n",
                        collapse("  real<lower=", lower_bound, "> ", 
                                 "prior_", pars[bound], "; \n"))
      out$model <- paste0("  // additionally draw samples from priors \n",
                          collapse("  prior_", pars[bound] ," ~ ",
                                   dis[bound], args[bound]," \n"))
    }
    if (any(!bound)) {  
      # unbounded parameters can be sampled in the generatated quantities block
      if (!is.null(hs_df)) {
        args[match("b", pars)] <- "(0, prior_hs_local * prior_hs_global);" 
      } 
      out$genD <- collapse("  real prior_", pars[!bound], "; \n")
      out$genC <- paste0("  // additionally draw samples from priors \n",
                         collapse("  prior_", pars[!bound], " <- ",
                                  dis[!bound], "_rng", args[!bound], " \n"))
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

stan_eta_re <- function(ranef, par = "") {
  # Write the random effects part of the linear predictor
  # Args:
  #   ranef: a named list returned by gather_ranef
  #   par: an optional character string to add to the variable names
  #        (used for non-linear models)
  # Returns:
  #   A string containing the random effects part of the linear predictor
  eta_re <- ""
  for (i in seq_along(ranef)) {
    pi <- if (nchar(par)) paste0(par, "_", i) else i
    if (length(ranef[[i]]) == 1 || attr(ranef[[i]], "cor")) {
      eta_re <- paste0(eta_re, " + Z_", pi,"[n] * r_", pi,"[J_", pi,"[n]]")
    } else {
      k <- seq_along(ranef[[i]])
      eta_re <- paste0(eta_re, collapse(" + Z_", pi, "[n, ", k, "]",
                                        " * r_", pi, "_" , k, "[J_", pi,"[n]]"))
    }
  }
  eta_re
}

stan_eta_transform <- function(family, link, add = FALSE) {
  # indicate whether eta needs to be transformed
  # Args:
  #   add: is the model weighted, censored, or truncated?
  !(add || family == "gaussian" && link == "log"
    || is.ordinal(family) || is.categorical(family) 
    || is.count(family) && link == "log" 
    || is.binary(family) && link == "logit"
    || is.zero_inflated(family) || is.hurdle(family))
}

stan_eta_ilink <- function(family, link) {
  # correctly apply inverse link to eta
  ilink <- stan_ilink(link)
  fl <- ifelse(family %in% c("gamma", "exponential"), 
               paste0(family,"_",link), family)
  switch(fl, c(paste0(ilink,"("), ")"),
         gamma_log = c("shape * exp(-(", "))"),
         gamma_inverse = c("shape * (", ")"),
         gamma_identity = c("shape / (", ")"),
         exponential_log = c("exp(-(", "))"),
         exponential_inverse = c("(", ")"),
         exponential_identity = c("inv(", ")"),
         weibull = c(paste0(ilink,"(("), ") / shape)"))
}
