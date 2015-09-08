# This file contains all functions generating Stan code

stan_model <- function(formula, data = NULL, family = "gaussian", link = "identity",
                       prior = prior_frame(), partial = NULL, threshold = "flexible", 
                       cov.ranef = NULL, sample.prior = FALSE, autocor = cor_arma(), 
                       save.model = NULL) {
  # Writes the regression model in Stan language
  # 
  # Args: 
  #   As in brm
  #
  # Returns:
  #  The model in stan code
  ee <- extract_effects(formula = formula, family = family, partial = partial) 
  if (family == "gaussian" && length(ee$response) > 1)
    family <- "multinormal"
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat") 
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  is_count <- family %in% c("poisson", "negbinomial", "geometric")
  is_multi <- family == "multinormal"

  if (family == "categorical") {
    X <- data.frame()
    Xp <- get_model_matrix(ee$fixed, data, rm_intercept = is_ordinal)
  } else {
    X <- get_model_matrix(ee$fixed, data, rm_intercept = is_ordinal)
    Xp <- get_model_matrix(partial, data, rm_intercept = TRUE)
  }  
  fixef <- colnames(X)
  paref <- colnames(Xp)
  Z <- lapply(ee$random, get_model_matrix, data = data)
  ranef <- lapply(Z,colnames)
  trait <- ifelse(is_multi, "_trait", "")
  names_cov_ranef <- names(cov.ranef)
  
  if (length(ee$group)) {
    # call stan_ranef for each random term seperately
    text_ranef <- lapply(1:length(ee$group), stan_ranef, 
                          ranef = ranef, group = ee$group, 
                          cor = ee$cor, prior = prior, 
                          names_cov_ranef = names_cov_ranef)
  } else {
    text_ranef <- list()
  }
  # combine random effects stan code of different grouping factors by names
  text_ranef <- collapse_lists(text_ranef)
  
  # generate important parts of the stan code
  text_eta <- stan_eta(family = family, link = link, fixef = fixef, 
                       paref = paref,  group = ee$group, autocor = autocor)
  text_ma <- stan_ma(family = family, link = link, autocor = autocor)
  text_ordinal <- stan_ordinal(family = family, link = link, 
                               partial = length(paref), 
                               threshold = threshold)  
  text_multi <- stan_multi(family, response = ee$response)
  text_llh <- stan_llh(family, link = link, 
                       add = is.formula(ee[c("se", "trials")]), 
                       weights = is.formula(ee$weights), 
                       cens = is.formula(ee$cens))
  if (is.formula(ee$cens) || is.formula(ee$weights) 
      || is_ordinal || family == "categorical") {
    text_llh <- paste0("  for (n in 1:N",trait,") { \n  ",text_llh,"  } \n")
  }
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    if (length(fixef)) 
      stan_prior(class = "b", coef = fixef, prior = prior),
    if (length(paref)) 
      stan_prior(class = "bp", coef = paref, prior = prior),
    if (autocor$p) 
      stan_prior(class = "ar", prior = prior),
    if (autocor$q) 
      stan_prior(class = "ma", prior = prior),
    if (is_ordinal && threshold == "flexible") 
      stan_prior("b_Intercept", prior = prior)
    else if (is_ordinal && threshold == "equidistant") 
      paste0(stan_prior(class = "b_Intercept1", prior = prior), 
             stan_prior(class = "delta", prior = prior)),
    if (family %in% c("gamma", "weibull", "negbinomial")) 
      stan_prior(class = "shape", prior = prior),
    if (family == "student") 
      stan_prior(class = "nu", prior = prior),
    if (is_linear && !is.formula(ee$se)) 
      stan_prior(class = "sigma", coef = ee$response, prior = prior), 
    if (is_multi) 
      paste0(stan_prior(class = "sigma", coef = ee$response, prior = prior),
             stan_prior(class = "Lrescor", prior = prior)),
    text_ranef$model)
  
  # generate code to additionally sample from priors if sample.prior = TRUE
  text_rngprior <- stan_rngprior(sample.prior = sample.prior, prior = text_prior, family = family)
  
  kronecker <- any(sapply(mapply(list, ranef, ee$group, SIMPLIFY = FALSE), 
                          function(x, names) length(x[[1]]) > 1 && x[[2]] %in% names, 
                          names = names_cov_ranef))
  text_functions <- stan_function(kronecker = kronecker)
  
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N; \n", 
    if (is_linear || is_skew) 
      "  real Y[N]; \n"
    else if (family %in% c("binomial", "bernoulli", "categorical")
             || is_count || is_ordinal) 
      "  int Y[N]; \n"
    else if (is_multi) 
      paste0("  int<lower=1> N_trait; \n  int<lower=1> K_trait; \n",  
             "  int NC_trait; \n  vector[K_trait] Y[N_trait]; \n"),
    if (length(fixef)) 
      "  int<lower=1> K; \n  matrix[N,K] X; \n",
    if (length(paref)) 
      "  int<lower=1> Kp; \n  matrix[N,Kp] Xp; \n",  
    if (autocor$p && is(autocor, "cor_arma")) 
      "  int<lower=1> Kar; \n  matrix[N,Kar] Yar; \n",
    if (autocor$q && is(autocor, "cor_arma")) 
      "  int<lower=1> Kma; \n  row_vector[Kma] Ema_pre[N]; \n  vector[N] tgroup; \n",
    if (is_linear && is.formula(ee$se))
      "  real<lower=0> sigma[N]; \n",
    if (is.formula(ee$weights))
      paste0("  vector<lower=0>[N",trait,"] weights; \n"),
    if (family == "binomial")
      paste0("  int trials", if (is.formula(ee$trials)) "[N]", "; \n"),
    if (is_ordinal || family == "categorical")
      paste0("  int ncat; \n"),
    if (is.formula(ee$cens) && !(is_ordinal || family == "categorical"))
      "  vector[N] cens; \n",
    text_ranef$data,
    "} \n")
  
  if (family == "categorical")
    zero <- c("  row_vector[1] zero; \n", "  zero[1] <- 0; \n")
  else zero <- NULL
  text_transformed_data <- paste0(
    "transformed data { \n",
      zero[1], text_ranef$tdataD, 
      zero[2], text_ranef$tdataC,
    "} \n")
  
  text_parameters <- paste0(
    "parameters { \n",
    if (length(fixef)) 
      "  vector[K] b; \n",
    if (length(paref)) 
      paste0("  matrix[Kp,ncat-1] bp; \n"),
    text_ordinal$par, text_ranef$par,
    if (autocor$p && is(autocor, "cor_arma")) 
      "  vector[Kar] ar; \n",
    if (autocor$q && is(autocor, "cor_arma")) 
      "  vector[Kma] ma; \n",
    if (is_linear && !is.formula(ee$se)) 
      "  real<lower=0> sigma; \n",
    if (family == "student") 
      "  real<lower=0> nu; \n",
    if (is_multi) 
      paste0("  vector<lower=0>[K_trait] sigma; \n",
             "  cholesky_factor_corr[K_trait] Lrescor; \n"),
    if (family %in% c("gamma", "weibull", "negbinomial")) 
      "  real<lower=0> shape; \n",
    text_rngprior$par,
    "} \n")
  
  # loop over all observations in transformed parameters if necessary
  make_loop <- length(ee$group) || autocor$q || text_eta$transform ||
                 (is_ordinal && !(family == "cumulative" && link == "logit"))
  if (make_loop && !is_multi) {
    text_loop <- c("  for (n in 1:N) { \n", "  } \n")
  } else if (is_multi) {
    text_loop <- c(paste0("  for (m in 1:N_trait) { \n  for (k in 1:K_trait) { \n",    
                           "    int n; \n    n <- (k-1)*N_trait + m; \n"), "  }} \n")
  } else {
    text_loop <- rep("", 2)
  }

  # combine all elements into on stan model
  model <- paste0(
    text_functions,
    text_data, 
    text_transformed_data, 
    text_parameters,
  "transformed parameters { \n",
    text_eta$transD, 
    text_ma$transD, 
    text_ordinal$transD, 
    text_ranef$transD, 
    text_eta$transC1, 
    text_ma$transC1, 
    text_ordinal$transC1, 
    text_ranef$transC, 
    text_loop[1],
      text_eta$transC2, 
      text_ma$transC2, 
      text_ordinal$transC2, 
      text_eta$transC3, 
    text_loop[2],
  "} \n",
  "model { \n",
    if (is.formula(ee$weights) && !is.formula(ee$cens)) 
      paste0("  vector[N",trait,"] lp_pre; \n"),
    text_prior, 
    text_llh, 
    if (is.formula(ee$weights) && !is.formula(ee$cens)) 
    "  increment_log_prob(dot_product(weights,lp_pre)); \n",
    text_rngprior$model,
  "} \n",
  "generated quantities { \n",
    text_multi$genD, 
    text_ranef$genD, 
    text_rngprior$genD, 
    text_multi$genC, 
    text_ranef$genC, 
    text_rngprior$genC,
  "} \n")
  
  # write the stan code to a file if save.model is a character string
  class(model) <- c("character", "brmsmodel")
  if (is.character(save.model)) {
    sink(save.model)
    cat(model)
    sink()
  }
  model
}

stan_ranef <- function(i, ranef, group, cor, prior = list(), 
                       names_cov_ranef = NULL) {
  # Random effects in Stan 
  # 
  # Args:
  #   i: the index of the grouping factor
  #   ranef: a list of random effects 
  #   group: a vector of grouping factors
  #   cor: a logical vector to indicate if correlations should be estimated
  #   prior: user defined priors
  #   names_cov_ranef: names of the grouping factors for which custom covariance matrices are specified.
  #
  # Returns:
  #   A vector of strings containing the random effects in stan language
  r <- ranef[[i]]
  g <- group[[i]]
  cor <- cor[[i]]
  ccov <- g %in% names_cov_ranef
  out <- list()
  out$data <- paste0("  int<lower=1> lev_",i,"[N]; \n",
                     "  int<lower=1> N_",i,"; \n",
                     "  int<lower=1> K_",i,"; \n",
                     if (ccov && (cor || length(r) == 1)) 
                       paste0("  matrix[N_",i,", N_",i,"] cov_",i,"; \n"),
                     if (ccov && !cor && length(r) > 1) 
                       paste0("  matrix[N_",i,"*K_",i,", N_",i,"*K_",i,"] cov_",i,"; \n"))
  out$model <- stan_prior(class = "sd", group = i, coef = r, prior = prior)
                      
  if (length(r) == 1) {  # only one random effect
    out$data <- paste0(out$data, "  real Z_",i,"[N]; \n")
    out$par <- paste0("  vector[N_",i,"] pre_",i,"; \n",
                      "  real<lower=0> sd_",i,"; \n")
    out$model <- paste0(out$model,"  pre_",i," ~ normal(0,1); \n")
    out$transD <- paste0("  vector[N_",i,"] r_",i,"; \n")
    out$transC <- paste0("  r_",i, " <- sd_",i," * (", 
                         if (ccov) paste0("cov_",i," * "), "pre_",i,"); \n")
  }  
  else if (length(r) > 1) {  # multiple random effects
    out$data <- paste0(out$data,  "  row_vector[K_",i,"] Z_",i,"[N]; \n  int NC_",i,"; \n")
    out$par <- paste0("  matrix[N_",i,",K_",i,"] pre_",i,"; \n",
                      "  vector<lower=0>[K_",i,"] sd_",i,"; \n",
                      if (cor) paste0("  cholesky_factor_corr[K_",i,"] L_",i,"; \n"))
    out$model <- paste0(out$model, 
                        if (cor) stan_prior(class = "L", group = i, prior = prior),
                        "  to_vector(pre_",i,") ~ normal(0,1); \n")
    out$transD <- paste0("  vector[K_",i,"] r_",i,"[N_",i,"]; \n")
    if (ccov) {  # customized covariance matrix supplied
      if (cor) {  # estimate correlations between random effects
        out$transC <- paste0("  r_",i," <- to_array(kronecker_cholesky(cov_",i,", L_",i,", sd_",i,") * ",
                             "to_vector(pre_",i,"), N_",i,", K_",i,"); \n")
      } else { 
        out$transC <- paste0("  r_",i," <- to_array(to_vector(rep_matrix(sd_",i,", N_",i,")) .* ",
                                "(cov_",i," * to_vector(pre_",i,")), N_",i,", K_",i,"); \n")
      }
    } else { 
      out$transC <- paste0("  for (i in 1:N_",i,") { \n",
                           "    r_",i, "[i] <- sd_",i," .* (", 
                           if (cor) paste0("L_",i," * "), 
                           "to_vector(pre_",i,"[i])); \n  } \n")
    }
    if (cor) {  # return correlations above the diagonal only
      out$genD <- paste0("  corr_matrix[K_",i,"] Cor_",i,"; \n",
                         "  vector<lower=-1,upper=1>[NC_",i,"] cor_",i,"; \n")
      out$genC <- paste0("  Cor_",i," <- multiply_lower_tri_self_transpose(L_",i,"); \n",
                         collapse(unlist(lapply(2:length(r), function(k) lapply(1:(k-1), function(j)
                           paste0("  cor_",i,"[",(k-1)*(k-2)/2+j,"] <- Cor_",i,"[",j,",",k,"]; \n")))))) 
    }  
  }
  out
}

stan_llh <- function(family, link, add = FALSE, 
                     weights = FALSE, cens = FALSE) {
  # Likelihoods in stan language
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   add: logical; indicating if there is information on se, trials or cat
  #   weights: logical; weights present?
  #   cens: logical; censored data?
  #
  # Returns:
  #   a string containing the likelihood of the model in stan language
  is_cat <- family %in% c("cumulative", "cratio", "sratio", "acat", "categorical")
  is_count <- family %in% c("poisson","negbinomial", "geometric")
  is_skew <- family %in% c("gamma","exponential","weibull")
  is_binary <- family %in% c("binomial", "bernoulli")
  
  simplify <- !cens && (is_binary && link == "logit" || is_count && link == "log" ||
                family %in% c("cumulative", "categorical") && link == "logit" && !add) 
  n <- ifelse(cens || weights || is_cat, "[n]", "")
  ns <- ifelse(add && (cens || weights), "[n]", "")
  ilink <- ifelse(cens && (is_binary && link == "logit" || is_count && link == "log"), 
                  stan_ilink(link), "")
  
  lin.args <- paste0("eta",n,",sigma",ns)
  if (simplify) { 
    llh.pre <- switch(family,
      poisson = c("poisson_log", paste0("eta",n)), 
      negbinomial = c("neg_binomial_2_log", paste0("eta",n,",shape")),
      geometric = c("neg_binomial_2_log", paste0("eta",n,",1")),
      cumulative = c("ordered_logistic", "eta[n],b_Intercept"),
      categorical = c("categorical_logit", 
                    "to_vector(append_col(zero, eta[n] + etap[n]))"), 
      binomial = c("binomial_logit", paste0("trials",ns,",eta",n)), 
      bernoulli = c("bernoulli_logit", paste0("eta",n)))
  } else {
    llh.pre <- switch(ifelse(is_cat, "categorical", 
                             ifelse(family == "gaussian" && link == "log", 
                                    "lognormal", family)),
      gaussian = c("normal", lin.args),
      student = c("student_t", paste0("nu,",lin.args)),
      cauchy = c("cauchy", lin.args),
      multinormal = c("multi_normal_cholesky", 
                      paste0("etam",n,",diag_pre_multiply(sigma,Lrescor)")),
      lognormal = c("lognormal", lin.args),
      poisson = c("poisson", paste0(ilink,"(eta",n,")")),
      negbinomial = c("neg_binomial_2", paste0(ilink,"(eta",n,"),shape")),
      geometric = c("neg_binomial_2", paste0(ilink,"(eta",n,"),1")),
      binomial = c("binomial", paste0("trials",ns,",",ilink,"(eta",n,")")),
      bernoulli = c("bernoulli", paste0(ilink,"(eta",n,")")), 
      gamma = c("gamma", paste0("shape,eta",n)), 
      exponential = c("exponential", paste0("eta",n)),
      weibull = c("weibull", paste0("shape,eta",n)), 
      categorical = c("categorical","p[n]"))
  }
  
  # write likelihood code
  type <- c("cens", "weights")[match(TRUE, c(cens, weights))]
  if (is.na(type)) type <- "general"
  addW <- ifelse(weights, "weights[n] * ", "")
  llh <- switch(type, 
    cens = paste0("if (cens[n] == 0) ", 
      ifelse(!weights, paste0("Y[n] ~ ", llh.pre[1],"(",llh.pre[2],"); \n"),
             paste0("increment_log_prob(", addW, llh.pre[1], "_log(Y[n],",llh.pre[2],")); \n")),
      "    else { \n",         
      "      if (cens[n] == 1) increment_log_prob(", addW, llh.pre[1], "_ccdf_log(Y[n],", llh.pre[2],")); \n",
      "      else increment_log_prob(", addW, llh.pre[1], "_cdf_log(Y[n],", llh.pre[2],")); \n",
      "    } \n"),
    weights = paste0("  lp_pre[n] <- ", llh.pre[1], "_log(Y[n],",llh.pre[2],"); \n"),
    general = paste0("  Y", n, " ~ ", llh.pre[1],"(",llh.pre[2],"); \n")) 
  llh
}

stan_eta <- function(family, link, fixef, paref = NULL, 
                     group = NULL, autocor = cor_arma()) {
  # linear predictor in Stan
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   fixef: names of the fixed effects parameters
  #   paref: names of the partiel effects paraneters
  #   group: names of the grouping factors
  #   autocor: autocorrelation structure
  # 
  # Return:
  #   the linear predictor in stan language
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat") 
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  is_count <- family %in% c("poisson", "negbinomial", "geometric")
  is_binary <- family %in% c("binomial", "bernoulli")
  is_multi <- family == "multinormal"
  
  eta <- list()
  # initialize eta
  eta$transD <- paste0("  vector[N] eta; \n", 
                       ifelse(length(paref), paste0("  matrix[N,ncat-1] etap; \n"), ""),
                       ifelse(is_multi, "  vector[K_trait] etam[N_trait]; \n", ""))
  eta.multi <- ifelse(is_multi, "etam[m,k]", "eta[n]")
  
  # transform eta before it is passed to the likelihood
  ilink <- stan_ilink(link)
  eta$transform <- !(link == "identity" || family == "gaussian" && link == "log" ||
                     is_ordinal || family == "categorical" || is_count && link == "log" ||
                     is_binary && link == "logit")
  eta_ilink <- rep("", 2)
  if (eta$transform) {
    eta_ilink <- switch(family, c(paste0(ilink,"("), ")"),
                   gamma = c(paste0("shape/(",ilink,"("), "))"), 
                   exponential = c(paste0(ilink,"(-("), "))"), 
                   weibull = c(paste0("inv(",ilink,"(-("), ")/shape))"))
    if (autocor$q > 0) {
      eta$transC3 <- paste0("    ",eta.multi," <- ",eta_ilink[1], eta.multi, eta_ilink[2],"; \n")
      eta_ilink <- rep("", 2)  
    }
  }
  
  # define fixed, random, and autocorrelation effects
  eta$transC1 <- paste0("  eta <- ", ifelse(length(fixef), "X*b", "rep_vector(0,N)"), 
                        if (autocor$p && is(autocor, "cor_arma")) " + Yar*ar", "; \n", 
                        if (length(paref)) "  etap <- Xp * bp; \n")
  if (length(group)) {
    ind <- 1:length(group)
    eta.re <- collapse(" + Z_",ind,"[n]*r_",ind,"[lev_",ind,"[n]]")
  } else {
    eta.re <- ""
  }
  eta.ma <- ifelse(autocor$q && is(autocor, "cor_arma"), " + Ema[n]*ma", "")
  if (nchar(eta.re) || nchar(eta.ma) || is_multi || nchar(eta_ilink[1])) {
    eta$transC2 <- paste0("    ",eta.multi," <- ",
                          eta_ilink[1],"eta[n]", eta.ma, eta.re, eta_ilink[2],"; \n")
  }
  eta
}

stan_ma <- function(family, link, autocor) {
  # moving average autocorrelation in Stan
  # 
  # Args:
  #   family: the model family
  #   link: the link function
  #   autocor: autocorrelation structure
  #
  # Returns:
  #   stan code for computing moving average effects
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_multi <- family == "multinormal"
  ma <- list()
  if (is(autocor, "cor_arma") && autocor$q) {
    link.fun <- c(identity = "", log = "log", inverse = "inv")[link]
    if (!(is_linear || is_multi))
      stop(paste("moving-average models for family", family, "are not yet implemented"))
    index <- ifelse(is_multi, "m,k", "n")
    ma$transD <- paste0("  row_vector[Kma] Ema[N]; \n  vector[N] e; \n") 
    ma$transC1 <- "  Ema <- Ema_pre; \n" 
    ma$transC2 <- paste0("    e[n] <- ",link.fun,"(Y[",index,"]) - eta[n]", "; \n", 
                         "    for (i in 1:Kma) if (n+1-i > 0 && n < N && tgroup[n+1] == tgroup[n+1-i]) \n",
                         "      Ema[n+1,i] <- e[n+1-i]", "; \n")
  }
  ma
}

stan_function <- function(kronecker = FALSE) {
  # stan code for user defined functions
  #
  # Args:
  #   kronecker: logical; is the kronecker product needed?
  #
  # Returns:
  #   a string containing defined functions in stan code
  out <- NULL
  if (kronecker) out <- paste0(
    "  // calculate the cholesky factor of the kronecker covariance matrix \n",
    "  matrix kronecker_cholesky(matrix X, matrix L, vector sd) { \n",
    "    matrix[rows(X)*rows(L), cols(X)*cols(L)] kron; \n",
    "    matrix[rows(L), cols(L)] C; \n",
    "    int rX; \n",
    "    int rC; \n",
    "    C <- multiply_lower_tri_self_transpose(L); \n",
    "    rX <- rows(X); \n",
    "    rC <- rows(C); \n",
    "    for (i in 1:rX) { \n",
    "      for (j in 1:rC) { \n",
    "        for (k in 1:rX) { \n",
    "          for (l in 1:rC) { \n",
    "            kron[(k-1)*rC+l, (i-1)*rC+j] <- sd[l] * sd[j] * X[k,i] * C[l,j]; \n",
    "          } \n",
    "        } \n",
    "      } \n",
    "    } \n",
    "    return cholesky_decompose(kron); \n",
    "  } \n",
    "  // turn a vector into a 2 dimensional array \n",
    "  vector[] to_array(vector X, int N, int K) { \n",
    "    vector[K] Y[N]; \n",
    "    for (i in 1:N) \n",
    "      Y[i] <- segment(X, (i-1)*K+1, K); \n",
    "    return Y; \n",
    "  } \n")
  return(paste0("functions { \n", out, "} \n"))
}

stan_multi <- function(family, response) {
  # multinormal effects in Stan
  #
  # Args:
  #   family: the model family
  #   response: names of the response variables
  # 
  # Returns: 
  #   string containing the stan code specific for multinormal models
  out <- list()
  if (family == "multinormal") {
   out$genD <- paste0("  matrix[K_trait,K_trait] Rescor; \n",
    "  vector<lower=-1,upper=1>[NC_trait] rescor; \n")
   out$genC <- paste0("  Rescor <- multiply_lower_tri_self_transpose(Lrescor); \n",
        collapse(unlist(lapply(2:length(response), function(i) lapply(1:(i-1), function(j)
        paste0("  rescor[",(i-1)*(i-2)/2+j,"] <- Rescor[",j,",",i,"]; \n"))))))
  }
  out
}

stan_ordinal <- function(family, link, partial = FALSE, threshold = "flexible") {
  # Ordinal effects in Stan
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   partial: logical; are there partial effects?
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  is_ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat")
  if (!(is_ordinal || family == "categorical")) return(list())
  ilink <- c(identity = "", log = "exp", inverse = "inv", 
             sqrt = "square", logit = "inv_logit", 
             probit = "Phi", probit_approx = "Phi_approx", 
             cloglog = "inv_cloglog")[link]
  th <- function(k) {
    sign <- ifelse(family %in% c("cumulative", "sratio")," - ", " + ")
    ptl <- ifelse(partial, paste0(sign, "etap[n,k]"), "") 
    if (sign == " - ") {
      out <- paste0("b_Intercept[",k,"]", ptl, " - eta[n]")
    } else {
      out <- paste0("eta[n]", ptl, " - b_Intercept[",k,"]")
    }
  }  
  sc <- ifelse(family == "sratio", "1-", "")
  intercept <- paste0("  ", ifelse(family == "cumulative", "ordered", "vector"), 
                      "[ncat-1] b_Intercept; \n")
  
  out <- list()
  if (is_ordinal) {
    if (threshold == "flexible") {
      out$par <- intercept
    } else if (threshold == "equidistant") {
      out$par <- paste0("  real b_Intercept1; \n",
                        "  real", if (family == "cumulative") "<lower=0>",
                        " delta; \n")
      out$transC1 <- paste0("  for (k in 1:(ncat-1)) { \n",
                            "    b_Intercept[k] <- b_Intercept1 + (k-1.0)*delta; \n",
                            "  } \n")
      out$transD <- intercept
    }
  }
  if (!(family %in% c("cumulative", "categorical") && ilink == "inv_logit")) {
    out$transD <- paste0(out$transD, "  vector[ncat] p[N]; \n", 
                         if (!family %in% c("cumulative", "categorical")) 
                           paste0("  vector[ncat-1] q[N]; \n"))
    if (family == "categorical" && ilink == "inv_logit") {
      out$transC <- paste0(
      "    p[n,1] <- 1.0; \n",
      "    for (k in 2:ncat) { \n",
      "      p[n,k] <- exp(eta[n,k-1]); \n",
      "    } \n",
      "    p[n] <- p[n]/sum(p[n]); \n")
    } else if (family == "cumulative") {
      out$transC2 <- paste0(
      "    p[n,1] <- ",ilink,"(",th(1),"); \n",
      "    for (k in 2:(ncat-1)) { \n", 
      "      p[n,k] <- ",ilink,"(",th("k"),") - ",ilink,"(",th("k-1"),"); \n", 
      "    } \n",
      "    p[n,ncat] <- 1 - ",ilink,"(",th("ncat-1"),"); \n")
    } else if (family %in% c("sratio", "cratio")) {
      out$transC2 <- paste0(
      "    for (k in 1:(ncat-1)) { \n",
      "      q[n,k] <- ",sc, ilink,"(",th("k"),"); \n",
      "      p[n,k] <- 1-q[n,k]; \n",
      "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n", 
      "    } \n",
      "    p[n,ncat] <- prod(q[n]); \n")
    } else if (family == "acat") {
      if (ilink == "inv_logit") {
        out$transC2 <- paste0(
        "    p[n,1] <- 1.0; \n",
        "    for (k in 1:(ncat-1)) { \n",
        "      q[n,k] <- ",th("k"),"; \n",
        "      p[n,k+1] <- q[n,1]; \n",
        "      for (kk in 2:k) p[n,k+1] <- p[n,k+1] + q[n,kk]; \n",
        "      p[n,k+1] <- exp(p[n,k+1]); \n",
        "    } \n",
        "    p[n] <- p[n]/sum(p[n]); \n")
      } else {
        out$transC2 <- paste0(                   
        "    for (k in 1:(ncat-1)) \n",
        "      q[n,k] <- ",ilink,"(",th("k"),"); \n",
        "    for (k in 1:ncat) { \n",     
        "      p[n,k] <- 1.0; \n",
        "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n",
        "      for (kk in k:(ncat-1)) p[n,k] <- p[n,k] * (1-q[n,kk]); \n",      
        "    } \n",
        "    p[n] <- p[n]/sum(p[n]); \n")
      }
    }
  }
  out
}

stan_prior <- function(class, coef = NULL, group = NULL, 
                       prior = prior_frame(), s = 2) {
  # Define priors for parameters in Stan language
  # 
  # Args:
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   prior: a data.frame containing user defined priors as returned by check_prior
  #   s: An integer >= 0 defining the number of spaces in front of the output string
  # 
  # Returns:
  #   A character strings in stan language that defines priors for a given class of parameters
  #   If a parameter has has no corresponding prior in prior 
  #   and also no internal default in stan_prior, an empty string is returned.
  
  # only consider user defined priors related to this class and group
  keep <- which(prior$class == class & (prior$coef %in% coef | !nchar(prior$coef)))
  user_prior <- prior[keep, ]
  if (!is.null(group)) {
    keep2 <- which(user_prior$group == group | !nchar(user_prior$group))
    user_prior <- user_prior[keep2, ]
  }
  
  # get base prior
  igroup <- which(with(user_prior, !nchar(coef) & nchar(group) & nchar(prior)))
  iclass <- which(with(user_prior, !nchar(coef) & !nchar(group) & nchar(prior)))
  if (length(igroup)) {  
    # if there is a global prior for this group
    base_prior <- user_prior[igroup, "prior"]
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
  s <- collapse(rep(" ", s))
  if (!is.null(group)) {
    class <- paste0(class,"_",group)
  }
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
  return(collapse(out))
}

stan_rngprior <- function(sample.prior, prior, family = "gaussian") {
  # stan code to sample from priors seperately
  #
  # Args:
  #   sample.prior: take samples from priors?
  #   prior: the character string taken from stan_prior
  #   family: the model family
  #
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  out <- list()
  if (sample.prior) {
    prior <- gsub(" ", "", paste0("\n",prior))
    pars <- gsub("\\\n|to_vector\\(|\\)", "", 
                 regmatches(prior, gregexpr("\\\n[^~]+", prior))[[1]])
    take <- !grepl("^pre_", pars)
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
    bound <- grepl("^sd|^sigma|^shape$|^nu$", pars) |  # do not change to ||
                   family == "cumulative" & grepl("^delta$", pars)
    if (any(bound)) {  
      # bounded parameters have to be sampled in the model block
      out$par <- collapse("  real<lower=0> prior_",pars[bound],"; \n")
      out$model <- collapse("  prior_",pars[bound]," ~ ",
                            dis[bound],args[bound]," \n")
    }
    if (any(!bound)) {  
      # unbounded parameters can be sampled in the generatated quantities block
      out$genD <- collapse("  real prior_",pars[!bound],"; \n")
      out$genC <- collapse("  prior_",pars[!bound]," <- ",
                           dis[!bound],"_rng",args[!bound]," \n")
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
         sqrt = "square", logit = "inv_logit", probit = "Phi", 
         probit_approx = "Phi_approx", cloglog = "inv_cloglog")
}
