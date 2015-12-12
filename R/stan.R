#' Stan Code for \pkg{brms} Models
#' 
#' Generate Stan code for \pkg{brms} models
#' 
#' @inheritParams brm
#' @param ... Other arguments for internal usage only
#' 
#' @return A character string containing the fully commented \pkg{Stan} code 
#'   to fit a \pkg{brms} model.
#'  
#' @examples 
#' make_stancode(rating ~ treat + period + carry + (1|subject), 
#'               data = inhaler, family = "cumulative")
#' 
#' make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c 
#'               + (1|patient) + (1|visit), 
#'               data = epilepsy, family = "poisson")
#'
#' @export
make_stancode <- function(formula, data = NULL, family = "gaussian", 
                          prior = NULL, autocor = NULL, 
                          multiply = NULL, partial = NULL, 
                          threshold = c("flexible", "equidistant"),
                          cov.ranef = NULL, sample.prior = FALSE, 
                          save.model = NULL, ...) {
  
  obj_family <- check_family(family) 
  link <- obj_family$link
  family <- obj_family$family
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  prior <- check_prior(prior, formula = formula, data = data, 
                       family = family, link = link, 
                       autocor = autocor, partial = partial, 
                       threshold = threshold) 
  et <- extract_time(autocor$formula)  
  ee <- extract_effects(formula, family = family, partial, et$all)
  data <- update_data(data, family = family, effects = ee, et$group)
  
  # flags to indicate of which type family is
  # see misc.R for details
  is_linear <- is.linear(family)
  is_ordinal <- is.ordinal(family)
  is_skewed <- is.skewed(family)
  is_count <- is.count(family)
  is_hurdle <- is.hurdle(family)
  is_zero_inflated <- is.zero_inflated(family)
  is_categorical <- family == "categorical"
  is_multi <- is_linear && length(ee$response) > 1
  has_sigma <- has_sigma(family, autocor = autocor, se = ee$se, 
                         is_multi = is_multi)
  has_shape <- has_shape(family)
  offset <- !is.null(model.offset(data)) 
  trunc <- get_boundaries(ee$trunc)  
  
  # generate fixed effects code
  if (is_categorical) {
    X <- data.frame()
    fixef <- colnames(X)
    Xp <- get_model_matrix(ee$fixed, data)
    temp_list <- check_intercept(colnames(Xp))
    paref <- temp_list$names
  } else {
    X <- get_model_matrix(ee$fixed, data)
    temp_list <- check_intercept(colnames(X))
    fixef <- temp_list$names
    Xp <- get_model_matrix(partial, data, rm_intercept = TRUE)
    paref <- colnames(Xp)
  }
  has_intercept <- temp_list$has_intercept
  multef <- colnames(get_model_matrix(multiply, data))
  text_fixef <- stan_fixef(fixef = fixef, multef = multef,
                           paref = paref, family = family, 
                           prior = prior, threshold = threshold,
                           has_intercept = has_intercept)
  
  # generate random effects code
  Z <- lapply(ee$random, get_model_matrix, data = data)
  ranef <- lapply(Z, colnames)
  trait <- ifelse(is_multi || is_hurdle || is_zero_inflated, "_trait", "")
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
  
  # generate other important parts of the stan code
  text_eta <- stan_eta(family = family, link = link, fixef = fixef, 
                       has_intercept = has_intercept, multef = multef,
                       paref = paref, group = ee$group, autocor = autocor,
                       add = is.formula(ee[c("weights", "cens", "trunc")]),
                       offset = offset, is_multi = is_multi)
  text_llh <- stan_llh(family, link = link, 
                       is_multi = is_multi,
                       se = is.formula(ee$se),  
                       weights = is.formula(ee$weights),
                       trials = is.formula(ee$trials),
                       cens = is.formula(ee$cens),
                       trunc = trunc, autocor = autocor)
  if (is.formula(ee$cens) || is.formula(ee$weights) || is.formula(ee$trunc) ||
      is_ordinal || is_categorical || is_hurdle || is_zero_inflated) {
    text_llh <- paste0("  for (n in 1:N",trait,") { \n  ",text_llh,"  } \n")
  }
  
  # generate stan code specific to certain models
  text_arma <- stan_arma(family = family, link = link, 
                         autocor = autocor, prior = prior,
                         se = is.formula(ee$se), is_multi = is_multi)
  text_multi <- stan_multi(family = family, response = ee$response,
                           prior = prior)
  text_ordinal <- stan_ordinal(family = family, link = link, 
                               prior = prior, partial = length(paref), 
                               threshold = threshold)  
  text_zi_hu <- stan_zero_inflated_hurdle(family = family)
  text_inv_gaussian <- stan_inv_gaussian(family = family, 
                                         weights = is.formula(ee$weights),
                                         cens = is.formula(ee$cens),
                                         trunc = is.formula(ee$trunc))
  kronecker <- needs_kronecker(names_ranef = ranef, names_group = ee$group,
                               names_cov_ranef = names_cov_ranef)
  text_misc_funs <- stan_misc_functions(link = link, kronecker = kronecker)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_fixef$prior,
    text_ordinal$prior,
    text_ranef$prior,
    text_arma$prior,
    text_multi$prior,
    if (has_sigma) 
      stan_prior(class = "sigma", coef = ee$response, prior = prior), 
    if (has_shape) 
      stan_prior(class = "shape", prior = prior),
    if (family == "student") 
      stan_prior(class = "nu", prior = prior),
    if (family == "beta") 
      stan_prior(class = "phi", prior = prior),
    stan_prior(class = "", prior = prior))
  # generate code to additionally sample from priors if sample.prior = TRUE
  text_rngprior <- stan_rngprior(sample.prior = sample.prior, 
                                 prior = text_prior, family = family,
                                 hs_df = attr(prior, "hs_df"))
  
  # generate functions block
  text_functions <- paste0(
    "functions { \n",
      text_arma$fun,
      text_zi_hu$fun,
      text_inv_gaussian$fun,
      text_misc_funs,
    "} \n")
  
  # generate data block
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  is_real_Y <- is_linear || is_skewed || 
               family %in% c("inverse.gaussian", "beta")
  is_int_Y <- family %in% c("binomial", "bernoulli", "categorical") || 
              is_count || is_ordinal
  N_bin <- ifelse(is.formula(ee$trials), "[N]", "")
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  # number of observations \n", 
    if (is_multi) {
      text_multi$data
    } else if (is_hurdle || is_zero_inflated) {
      text_zi_hu$data
    } else if (is_real_Y) {
      "  vector[N] Y;  # response variable \n"
    } else if (is_int_Y) {
      "  int Y[N];  # response variable \n"
    },
    text_fixef$data,
    text_ranef$data,
    text_arma$data,
    text_inv_gaussian$data,
    if (family %in% c("binomial", "zero_inflated_binomial"))
      paste0("  int trials", N_bin, ";  # number of trials \n"),
    if (is_ordinal || is_categorical)
      paste0("  int ncat;  # number of categories \n"),
    if (offset)
      "  vector[N] offset;  # added to the linear predictor \n",
    if (is.formula(ee$se) && !(use_cov(autocor) && (Kar || Kma)))
      "  vector<lower=0>[N] se;  # SEs for meta-analysis \n",
    if (is.formula(ee$weights))
      paste0("  vector<lower=0>[N",trait,"] weights;  # model weights \n"),
    if (is.formula(ee$cens))
      paste0("  vector[N",trait,"] cens;  # indicates censoring \n"),
    if (trunc$lb > -Inf)
      paste0("  ", ifelse(is_int_Y, "int", "real"), " lb;",  
             "  # lower bound for truncation; \n"),
    if (trunc$ub < Inf)
      paste0("  ", ifelse(is_int_Y, "int", "real"), " ub;",  
             "  # upper bound for truncation; \n"),
    "} \n")
  
  # generate transformed parameters block
  zero <- list()
  if (is_categorical) {
    zero$tdataD <- "  row_vector[1] zero; \n"
    zero$tdataC <- "  zero[1] <- 0; \n"
  }
  text_transformed_data <- paste0(
    "transformed data { \n",
       zero$tdataD,
       text_ranef$tdataD, 
       zero$tdataC,
       text_ranef$tdataC,
    "} \n")
  
  # generate parameters block
  text_parameters <- paste0(
    "parameters { \n",
    text_fixef$par,
    text_ordinal$par,
    text_ranef$par,
    text_arma$par,
    text_multi$par,
    if (has_sigma)
      "  real<lower=0> sigma;  # residual SD \n",
    if (family == "student") 
      "  real<lower=1> nu;  # degrees of freedom \n",
    if (has_shape) 
      "  real<lower=0> shape;  # shape parameter \n",
    if (family == "beta") 
      "  real<lower=0> phi;  # precision parameter \n",
    if (!is.null(attr(prior, "hs_df"))) 
      paste0("  # horseshoe shrinkage parameters \n",
             "  vector<lower=0>[K] hs_local; \n",
             "  real<lower=0> hs_global; \n"),
    text_rngprior$par,
    "} \n")
  
  # generate transformed parameters block
  # loop over all observations in transformed parameters if necessary
  make_loop <- length(ee$group) || (Kar || Kma) && !is.formula(ee$se) ||  
               text_eta$transform || 
               (is_ordinal && !(family == "cumulative" && link == "logit"))
  if (make_loop && !is_multi) {
    text_loop <- c(paste0("  # if available add REs to linear predictor \n",
                          "  for (n in 1:N) { \n"), "  } \n")
  } else if (is_multi) {
    text_loop <- text_multi$loop
  } else {
    text_loop <- rep("", 2)
  }
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_eta$transD, 
      text_arma$transD, 
      text_ordinal$transD,
      text_multi$transD,
      text_ranef$transD, 
      text_eta$transC1, 
      text_arma$transC1, 
      text_ordinal$transC1, 
      text_ranef$transC, 
      text_loop[1],
        text_eta$transC2, 
        text_arma$transC2, 
        text_ordinal$transC2, 
        text_eta$transC3, 
      text_loop[2],
      text_multi$transC,
    "} \n")
  
  # generate model block
  lp_pre_needed <- is.formula(ee$weights) && !is.formula(ee$cens)
  text_model <- paste0(
    "model { \n",
      if (lp_pre_needed) 
        paste0("  vector[N",trait,"] lp_pre; \n"),
      "  # prior specifications \n", 
      text_prior, 
      "  # likelihood contribution \n",
      text_llh, 
      if (lp_pre_needed)  
        "  increment_log_prob(dot_product(weights, lp_pre)); \n",
      text_rngprior$model,
    "} \n")
  
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_multi$genD, 
      text_ranef$genD, 
      text_rngprior$genD, 
      text_multi$genC, 
      text_ranef$genC, 
      text_rngprior$genC,
    "} \n")

  # combine all elements into a complete Stan model
  complete_model <- paste0(
    text_functions,
    text_data, 
    text_transformed_data, 
    text_parameters,
    text_transformed_parameters,
    text_model,
    text_generated_quantities)
  
  # write the stan code to a file if save.model is a character string
  class(complete_model) <- c("character", "brmsmodel")
  if (is.character(save.model)) {
    sink(save.model)
    cat(complete_model)
    sink()
  }
  complete_model
}

stan_fixef <- function(fixef, multef, paref, family = "gaussian", 
                       prior = prior_frame(), has_intercept = TRUE, 
                       threshold = "flexible") {
  # Stan code for fixec effects
  #
  # Args:
  #   fixef: names of the fixed effects
  #   multef: names of the multiplicative effects
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
  if (has_intercept && !is.ordinal(family)) {
    # intercepts for ordinal models are defined in stan_ordinal
    if (family == "categorical") {
      out$par <- paste0(out$par,
        "  row_vector[ncat - 1] b_Intercept;  # fixed effects intercepts \n")
    } else {
      out$par <- paste0(out$par,
        "  real b_Intercept;  # fixed effects Intercept \n")
    }
    out$prior <- paste0(out$prior, stan_prior("b_Intercept", prior = prior))
  }
  if (length(fixef)) {
    out$data <- paste0(out$data, 
      "  int<lower=1> K;  # number of fixed effects \n", 
      "  matrix[N, K] X;  # FE design matrix \n")
    out$par <- paste0(out$par,
      "  vector[K] b;  # fixed effects \n") 
    fixef_prior <- stan_prior(class = "b", coef = fixef, prior = prior)
    out$prior <- paste0(out$prior, fixef_prior)
  }
  if (length(multef)) {
    out$data <- paste0(out$data,
      "  int<lower=1> Km; \n",
      "  matrix[N, Km] Xm; \n")
    out$par <- paste0(out$par, 
      "  vector[Km] bm; \n")
    multef_prior <- stan_prior(class = "bm", coef = multef, prior = prior)
    out$prior <- paste0(out$prior, multef_prior)
  }
  if (length(paref)) {
    out$data <- paste0(out$data, 
     "  int<lower=1> Kp;  # number of category specific effects \n",
     "  matrix[N, Kp] Xp;  # CSE design matrix \n")
    out$par <- paste0(out$par,
     "  matrix[Kp, ncat - 1] bp;  # category specific effects \n")
    paref_prior <- stan_prior(class = "bp", coef = paref, prior = prior)
    out$prior <- paste0(out$prior, paref_prior)
  }
  out
}

stan_ranef <- function(i, ranef, group, cor, prior = prior_frame(), 
                       names_cov_ranef = NULL) {
  # Random effects in Stan 
  # 
  # Args:
  #   i: the index of the grouping factor
  #   ranef: a list of random effects 
  #   group: a vector of grouping factors
  #   cor: a logical vector to indicate if correlations should be estimated
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   names_cov_ranef: names of the grouping factors 
  #                    for which custom covariance matrices are specified.
  #
  # Returns:
  #   A vector of strings containing the random effects in stan language
  r <- ranef[[i]]
  g <- group[[i]]
  cor <- cor[[i]]
  ccov <- g %in% names_cov_ranef
  out <- list()
  out$data <- paste0("  # data for random effects of ",g," \n",
                     "  int<lower=1> J_",i,"[N];  # RE levels \n",
                     "  int<lower=1> N_",i,";  # number of levels \n",
                     "  int<lower=1> K_",i,";  # number of REs \n",
                     if (ccov && (cor || length(r) == 1)) 
                       paste0("  matrix[N_",i,", N_",i,"] cov_",i,";",
                              "  # user defined covariance matrix \n"),
                     if (ccov && !cor && length(r) > 1) 
                       paste0("  matrix[N_",i," * K_",i,", N_",i," * K_",i,"] cov_",i,";",
                              "  # user defined covariance matrix \n"))
  out$prior <- stan_prior(class = "sd", group = i, coef = r, prior = prior)
                      
  if (length(r) == 1) {  # only one random effect
    out$data <- paste0(out$data, "  real Z_",i,"[N];  # RE design matrix \n")
    out$par <- paste0("  vector[N_",i,"] pre_",i,";  # unscaled REs \n",
                      "  real<lower=0> sd_",i,";  # RE standard deviation \n")
    out$prior <- paste0(out$prior,"  pre_",i," ~ normal(0, 1); \n")
    out$transD <- paste0("  vector[N_",i,"] r_",i,";  # REs \n")
    out$transC <- paste0("  r_",i, " <- sd_",i," * (", 
                         if (ccov) paste0("cov_",i," * "), "pre_",i,");",
                         "  # scale REs \n")
  }  
  else if (length(r) > 1) {  # multiple random effects
    out$data <- paste0(out$data,  
      "  row_vector[K_",i,"] Z_",i,"[N];  # RE design matrix \n",  
      "  int NC_",i,";  # number of correlations \n")
    out$par <- paste0("  matrix[N_",i,", K_",i,"] pre_",i,";  # unscaled REs \n",
      "  vector<lower=0>[K_",i,"] sd_",i,";  # RE standard deviation \n",
      if (cor) paste0("  cholesky_factor_corr[K_",i,"] L_",i,
      ";  # cholesky factor of correlations matrix \n"))
    out$prior <- paste0(out$prior, 
      if (cor) stan_prior(class = "L", group = i, prior = prior),
      "  to_vector(pre_",i,") ~ normal(0, 1); \n")
    out$transD <- paste0("  vector[K_",i,"] r_",i,"[N_",i,"];  # REs \n")
    if (ccov) {  # customized covariance matrix supplied
      if (cor) {  # estimate correlations between random effects
        out$transC <- paste0(
          "  r_",i," <- to_array(kronecker_cholesky(cov_",i,", L_",i,", sd_",i,") * ",
          "to_vector(pre_",i,"), N_",i,", K_",i,");  # scale REs \n")
      } else { 
        out$transC <- paste0(
          "  r_",i," <- to_array(to_vector(rep_matrix(sd_",i,", N_",i,")) .* ",
          "(cov_",i," * to_vector(pre_",i,")), N_",i,", K_",i,");  # scale REs \n")
      }
    } else { 
      out$transC <- paste0(
        "  for (i in 1:N_",i,") { \n",
        "    r_",i, "[i] <- sd_",i," .* (", 
        if (cor) paste0("L_",i," * "), 
        "to_vector(pre_",i,"[i]));  # scale REs \n  } \n")
    }
    if (cor) {  # return correlations above the diagonal only
      cors_genC <- ulapply(2:length(r), function(k) lapply(1:(k-1), function(j)
        paste0("  cor_",i,"[",(k-1)*(k-2)/2+j,"] <- Cor_",i,"[",j,",",k,"]; \n")))
      out$genD <- paste0("  corr_matrix[K_",i,"] Cor_",i,"; \n",
                         "  vector<lower=-1,upper=1>[NC_",i,"] cor_",i,"; \n")
      out$genC <- paste0("  # take only relevant parts of correlation matrix \n",
                         "  Cor_",i," <- multiply_lower_tri_self_transpose(L_",i,"); \n",
                         collapse(cors_genC)) 
    }  
  }
  out
}

stan_llh <- function(family, link, se = FALSE, weights = FALSE, 
                     trials = FALSE, cens = FALSE, trunc = .trunc(), 
                     autocor = cor_arma(), is_multi = FALSE) {
  # Likelihoods in stan language
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   se: logical; user defined SEs present?
  #   weights: logical; weights present?
  #   trials: logical; number of bernoulli trials given per observation?
  #   cens: logical; censored data?
  #   trunc: list containing lower and upper truncation boundaries
  #   autocor: autocorrelation structure; an object of classe cor_arma
  #   multi: is the model multivariate?
  #
  # Returns:
  #   a string containing the likelihood of the model in stan language
  is_linear <- is.linear(family)
  is_catordinal <- is.ordinal(family) || family == "categorical"
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
      stop("Invalid addition arguments")
    }
  } else if (family == "gaussian" && link == "log") {
    # prepare for use of lognormal likelihood
    family <- "lognormal"
    link <- "identity"
  } 
  
  simplify <- !is_trunc && !cens && 
    (is_binary && link == "logit" || is_count && link == "log" ||
    family %in% c("cumulative", "categorical") && link == "logit") 
  n <- ifelse(cens || weights || is_trunc || is_catordinal ||
              is_hurdle || is_zero_inflated, "[n]", "")
  ns <- ifelse((se || trials) && (cens || weights || is_trunc) 
               || trials && is_zero_inflated, "[n]", "")
  sigma <- paste0(ifelse(se, "se", "sigma"), ns)
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
                    gamma_log = "shape * exp(-eta[n])",
                    gamma_inverse = "shape * eta[n]",
                    gamma_identity = "shape / eta[n]",
                    exponential_log = "exp(-eta[n])",
                    exponential_inverse = "eta[n]",
                    exponential_identity = "inv(eta[n])",
                    weibull = paste0(ilink,"(eta[n] / shape)"))
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
      cumulative = c("ordered_logistic", "eta[n], b_Intercept"),
      categorical = c("categorical_logit", 
                      "to_vector(append_col(zero, eta[n] + etap[n]))"), 
      binomial = c("binomial_logit", paste0("trials",ns,", eta",n)), 
      bernoulli = c("bernoulli_logit", paste0("eta",n)))
  } else {
    llh_pre <- switch(ifelse(is_catordinal, "categorical", family),
      gaussian = c("normal", paste0(eta,", ",sigma)),
      gaussian_cov = c("normal_cov", paste0(eta,", squared_se, N_tg, ", 
                       "begin_tg, nrows_tg, res_cov_matrix")),
      student = c("student_t",  paste0("nu, ",eta,", ",sigma)),
      student_cov = c("student_t_cov", paste0("nu, ",eta,", squared_se, N_tg, ", 
                      "begin_tg, nrows_tg, res_cov_matrix")),
      cauchy = c("cauchy", paste0(eta,", ", sigma)),
      cauchy_cov = c("student_t_cov", paste0("1, ",eta,", squared_se, N_tg, ", 
                     "begin_tg, nrows_tg, res_cov_matrix")),
      lognormal = c("lognormal", paste0(eta,", sigma",ns)),
      multi_gaussian = c("multi_normal_cholesky", paste0("Eta",n,", LSigma")),
      multi_student = c("multi_student_t", paste0("nu, Eta",n,", Sigma")),
      multi_cauchy = c("multi_student_t", paste0("1.0, Eta",n,", Sigma")),
      poisson = c("poisson", eta),
      negbinomial = c("neg_binomial_2", paste0(eta,", shape")),
      geometric = c("neg_binomial_2", paste0(eta,", 1")),
      binomial = c("binomial", paste0("trials",ns,", ",eta)),
      bernoulli = c("bernoulli", eta), 
      gamma = c("gamma", paste0("shape, ",eta)), 
      exponential = c("exponential", eta),
      weibull = c("weibull", paste0("shape, ",eta)), 
      inverse.gaussian = c("inv_gaussian", 
                           paste0(eta, ", shape, log_Y",n,", sqrt_Y",n)),
      beta = c("beta", paste0(eta, " * phi, (1 - ", eta, ") * phi")),
      categorical = c("categorical", "p[n]"),
      hurdle_poisson = c("hurdle_poisson", "eta[n], eta[n + N_trait]"),
      hurdle_negbinomial = c("hurdle_neg_binomial_2", 
                             "eta[n], eta[n + N_trait], shape"),
      hurdle_gamma = c("hurdle_gamma", "shape, eta[n], eta[n + N_trait]"),
      zero_inflated_poisson = c("zero_inflated_poisson", 
                                "eta[n], eta[n + N_trait]"),
      zero_inflated_negbinomial = c("zero_inflated_neg_binomial_2", 
                                    "eta[n], eta[n + N_trait], shape"),
      zero_inflated_binomial = c("zero_inflated_binomial", 
        paste0("trials",ns,", eta[n], eta[n + N_trait]")))
  }
  
  # write likelihood code
  type <- c("cens", "weights")[match(TRUE, c(cens, weights))]
  if (is.na(type)) type <- "general"
  # prepare for possible truncation
  code_trunc <- ""
  if (is_trunc) {
    if (type %in% c("cens", "weights")) {
      stop(paste("truncation is not yet possible in censored or weighted models"))
    } else {
      lb <- ifelse(trunc$lb > -Inf, "lb", "")
      ub <- ifelse(trunc$ub < Inf, "ub", "")
      code_trunc <- paste0(" T[",lb,", ",ub,"]")
    }
  }
  add_weights <- ifelse(weights, "weights[n] * ", "")
  llh <- switch(type, 
    cens = paste0("  # special treatment of censored data \n",
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
    weights = paste0("  lp_pre[n] <- ", llh_pre[1], "_log(Y[n], ",llh_pre[2],"); \n"),
    general = paste0("  Y", n, " ~ ", llh_pre[1],"(",llh_pre[2],")", 
                     code_trunc, "; \n")) 
  llh
}

stan_eta <- function(family, link, fixef, has_intercept = TRUE, 
                     multef = NULL, paref = NULL, group = NULL, 
                     autocor = cor_arma(),  add = FALSE, 
                     offset = FALSE, is_multi = FALSE) {
  # linear predictor in Stan
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   fixef: names of the fixed effects
  #   multef: names of the multiplicative effects 
  #   paref: names of the category specific effects
  #   group: names of the grouping factors
  #   autocor: autocorrelation structure
  #   add: is the model weighted, censored, or truncated?
  #   offset: is an offset defined?
  #   is_multi: is the model multivariate?
  # 
  # Return:
  #   the linear predictor in stan language
  is_linear <- is.linear(family)
  is_ordinal <- is.ordinal(family)
  is_cat <- family == "categorical"
  is_skewed <- is.skewed(family)
  is_count <- is.count(family) || is.zero_inflated(family) ||
              family %in% c("hurdle_poisson", "hurdle_negbinomial")
  is_binary <- is.binary(family) || family == "zero_inflated_binomial"
  
  eta <- list()
  # initialize eta
  eta$transD <- paste0(
    "  vector[N] eta;  # linear predictor \n", 
    if (length(multef)) 
      "  vector[N] etam;  # multiplicative linear predictor \n",
    if (length(paref)) 
      paste0("  matrix[N, ncat - 1] etap;",
             "  # linear predictor for category specific effects \n"),
    if (is_multi) 
      paste0("  vector[K_trait] Eta[N_trait];",
             "  # multivariate linear predictor matrix \n"))
  eta_obj <- ifelse(is_multi, "Eta[m, k]", "eta[n]")
  s <- ifelse(is_multi, "  ", "")
  
  # transform eta before it is passed to the likelihood
  ilink <- stan_ilink(link)
  eta$transform <- !(add || link == "identity"
                     || family == "gaussian" && link == "log"
                     || is_ordinal || family == "categorical" 
                     || is_count && link == "log" 
                     || is_binary && link == "logit"
                     || family == "hurdle_gamma")
  eta_ilink <- etam <- rep("", 2)
  if (length(multef)) etam <- c("etam[n] * (", ")")
  if (eta$transform || (get_ar(autocor) && !use_cov(autocor))) {
    fl <- ifelse(family %in% c("gamma", "exponential"), 
                 paste0(family,"_",link), family)
    eta_ilink <- switch(fl, c(paste0(ilink,"("), ")"),
                        gamma_log = c("shape * exp(-(", "))"),
                        gamma_inverse = c("shape * (", ")"),
                        gamma_identity = c("shape / (", ")"),
                        exponential_log = c("exp(-(", "))"),
                        exponential_inverse = c("(", ")"),
                        exponential_identity = c("inv(", ")"),
                        weibull = c(paste0(ilink,"(("), ") / shape)"))
    if (get_ar(autocor)) {
      eta_ar <- ifelse(!use_cov(autocor), " + head(E[n], Kar) * ar", "")
      eta$transC3 <- paste0("    ", s, eta_obj," <- ", eta_ilink[1], 
                            etam[1], eta_obj, eta_ar, etam[2], 
                            eta_ilink[2], "; \n")
      # don't apply link and multiplicative effects twice
      eta_ilink <- etam <- rep("", 2)
    }
  }
  
  # define fixed, random, and autocorrelation effects
  etap <- if (length(paref) || is_cat) {
    paste0("  etap <- ", 
           ifelse(length(paref), "Xp * bp", "rep_matrix(0, Kp, ncat - 1)"),
           if (is_cat && has_intercept) " + rep_matrix(b_Intercept, N)", "; \n")
  }
  if (length(group)) {
    ind <- 1:length(group)
    eta_re <- collapse(" + Z_",ind,"[n] * r_",ind,"[J_",ind,"[n]]")
  } else {
    eta_re <- ""
  }
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor), 
                   " + head(E[n], Kma) * ma", "")
  if (nchar(eta_re) || nchar(eta_ma) || is_multi || nchar(eta_ilink[1])) {
    eta$transC2 <- paste0("    ",s, eta_obj," <- ", eta_ilink[1], 
                          etam[1], "eta[n]", eta_ma, eta_re, etam[2],
                          eta_ilink[2],"; \n")
    eta_ilink <- etam <- rep("", 2)
  }
  # compute transC1 last to correctly incorporate multiplicative effects 
  eta$transC1 <- paste0(
    "  # compute linear predictor \n",
    if (length(multef)) "  etam <- exp(Xm * bm); \n",
    "  eta <- ", etam[1],
    ifelse(length(fixef), "X * b", "rep_vector(0, N)"), 
    if (has_intercept && !(is_ordinal || is_cat)) 
      " + b_Intercept",
    if (offset) " + offset",
    if (get_arr(autocor)) " + Yarr * arr", 
    etam[2], "; \n", etap)
  eta
}

stan_arma <- function(family, link, autocor, prior = prior_frame(),
                      se = FALSE, is_multi = FALSE) {
  # AR(R)MA autocorrelation in Stan
  # 
  # Args:
  #   family: the model family
  #   link: the link function
  #   autocor: autocorrelation structure; object of class cor_arma
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   se: user defined standard errors present?
  #   is_multi: is the model multivariate?
  #
  # Returns:
  #   stan code for computing AR(R)MA effects
  is_linear <- is.linear(family)
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  Karr <- get_arr(autocor)
  out <- list()
  if (Kar || Kma) {
    if (!is_linear) {
      stop(paste("ARMA effects for family", family, "are not yet implemented"))
    }
    out$data <- paste0(out$data,
      "  # data needed for ARMA effects \n",
      "  int<lower=0> Kar;  # AR order \n",
      "  int<lower=0> Kma;  # MA order \n",
      "  int<lower=1> Karma;  # max(Kma, Kar) \n",
      "  matrix[N, Karma] E_pre; # matrix of zeros \n",
      "  vector[N] tgroup;  # indicates independent groups \n")
    # restrict ARMA effects to be in [-1,1] when using covariance
    # formulation as they cannot be outside this interval anyway
    restrict <- ifelse(use_cov(autocor), "<lower=-1, upper=1>", "")
    if (Kar) {
      out$par <- paste0(out$par, 
        "  vector", restrict, "[Kar] ar;  # autoregressive effects \n")
      out$prior <- paste0(out$prior, stan_prior(class = "ar", prior = prior))
    }
    if (Kma) {
      out$par <- paste0(out$par, 
        "  vector", restrict, "[Kma] ma;  # moving-average effects \n")
      out$prior <- paste0(out$prior, stan_prior(class = "ma", prior = prior))
    }
    
    if (use_cov(autocor)) {
      # if the user wants ARMA effects to be estimated using
      # a covariance matrix for residuals
      if (is_multi) {
        stop(paste("multivariate models are not yet allowed", 
                   "when using ARMA covariance matrices"))
      }
      out$data <- paste0(out$data,
        "  # see the functions block for details \n",
        "  int<lower=1> N_tg; \n",   
        "  int<lower=1> begin_tg[N_tg]; \n",
        "  int<lower=1> nrows_tg[N_tg]; \n",
        "  vector[N] squared_se; \n")
      out$transD <- "  matrix[max(nrows_tg), max(nrows_tg)] res_cov_matrix; \n"
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
      out$transC1 <- paste0("  # compute residual covariance matrix; \n",
                            "  res_cov_matrix <- cov_matrix_", cov_mat_fun, 
                            "(", cov_mat_args, ", sigma, max(nrows_tg)); \n")
      # defined selfmade functions for the functions block
      if (family == "gaussian") {
        out$fun <- paste0(out$fun,
        "  /* multi normal log-PDF for special residual covariance structures \n",
        "   * currently only ARMA effects of order 1 are implemented \n",
        "   * Args: \n",
        "   *   y: response vector \n",
        "   *   eta: linear predictor \n",
        "   *   squared_se: square of the user defined standard errors \n",
        "   *               will be set to zero if non are defined \n",
        "   *   N_tg: number of groups \n",
        "   *   begin: indicates the first observation in each group \n",
        "   *   nrows: number of observations in each group \n",
        "   *   res_cov_matrix: AR1, MA1, or ARMA1 covariance matrix; \n",
        "   * Returns: \n",
        "   *   sum of the log-PDF values of all observations \n",
        "   */ \n",
        "   real normal_cov_log(vector y, vector eta, vector squared_se, \n", 
        "                        int N_tg, int[] begin, int[] nrows, \n",
        "                        matrix res_cov_matrix) { \n",
        "     vector[N_tg] log_post; \n",
        "     for (i in 1:N_tg) { \n",
        "       matrix[nrows[i], nrows[i]] Sigma; \n",
        "       vector[nrows[i]] y_part; \n",
        "       vector[nrows[i]] eta_part; \n",
        "       vector[nrows[i]] squared_se_part; \n",
        "       y_part <- segment(y, begin[i], nrows[i]); \n",
        "       eta_part <- segment(eta, begin[i], nrows[i]); \n",
        "       squared_se_part <- segment(squared_se, begin[i], nrows[i]); \n",
        "       Sigma <- block(res_cov_matrix, 1, 1, nrows[i], nrows[i]) \n",
        "                + diag_matrix(squared_se_part); \n",
        "       Sigma <- cholesky_decompose(Sigma); \n",
        "       log_post[i] <- multi_normal_cholesky_log(y_part, eta_part, Sigma); \n",
        "     } \n",                       
        "     return sum(log_post); \n",
        "   } \n")
      } else { # family %in% c("student", "cauchy")
        out$fun <- paste0(out$fun,
        "  /* multi student-t log-PDF for special residual covariance structures \n",
        "   * currently only ARMA effects of order 1 are implemented \n",
        "   * Args: \n",
        "   *   y: response vector \n",
        "   *   nu: degrees of freedom parameter \n",
        "   *   eta: linear predictor \n",
        "   *   squared_se: square of the user defined standard errors \n",
        "   *               will be set to zero if non are defined \n",
        "   *   N_tg: number of groups \n",
        "   *   begin: indicates the first observation in each group \n",
        "   *   nrows: number of observations in each group \n",
        "   *   res_cov_matrix: AR1, MA1, or ARMA1 covariance matrix; \n",
        "   * Returns: \n",
        "   *   sum of the log-PDF values of all observations \n",
        "   */ \n",
        "   real student_t_cov_log(vector y, real nu, vector eta, \n", 
        "                          vector squared_se, int N_tg, int[] begin, \n",
        "                          int[] nrows, matrix res_cov_matrix) { \n",
        "     vector[N_tg] log_post; \n",
        "     for (i in 1:N_tg) { \n",
        "       matrix[nrows[i], nrows[i]] Sigma; \n",
        "       vector[nrows[i]] y_part; \n",
        "       vector[nrows[i]] eta_part; \n",
        "       vector[nrows[i]] squared_se_part; \n",
        "       y_part <- segment(y, begin[i], nrows[i]); \n",
        "       eta_part <- segment(eta, begin[i], nrows[i]); \n",
        "       squared_se_part <- segment(squared_se, begin[i], nrows[i]); \n",
        "       Sigma <- block(res_cov_matrix, 1, 1, nrows[i], nrows[i]) \n",
        "                + diag_matrix(squared_se_part); \n",
        "       log_post[i] <- multi_student_t_log(y_part, nu, eta_part, Sigma); \n",
        "     } \n",                       
        "     return sum(log_post); \n",
        "   } \n")
      }
      if (Kar && !Kma) {
        out$fun <- paste0(out$fun,
        "  /* compute the covariance matrix for an AR1 process \n",
        "   * Args: \n",
        "   *   ar: AR1 autocorrelation \n",
        "   *   sigma: standard deviation of the AR1 process \n",
        "   *   nrows: number of rows of the covariance matrix \n",
        "   * Returns: \n",
        "   *   A nrows x nrows AR1 covariance matrix \n",
        "   */ \n",
        "   matrix cov_matrix_ar1(real ar, real sigma, int nrows) { \n",
        "     matrix[nrows, nrows] mat; \n",
        "     vector[nrows - 1] gamma; \n",
        "     mat <- diag_matrix(rep_vector(1, nrows)); \n",
        "     for (i in 2:nrows) { \n",
        "       gamma[i - 1] <- pow(ar, i - 1); \n",
        "       for (j in 1:(i - 1)) { \n",
        "         mat[i, j] <- gamma[i - j]; \n",
        "         mat[j, i] <- gamma[i - j]; \n",
        "       } \n",
        "     } \n",
        "     return sigma^2 / (1 - ar^2) * mat; \n",
        "   } \n")
      } else if (!Kar && Kma) {
        out$fun <- paste0(out$fun,
        "  /* compute the covariance matrix for an MA1 process \n",
        "   * Args: \n",
        "   *   ma: MA1 autocorrelation \n",
        "   *   sigma: standard deviation of the MA1 process \n",
        "   *   nrows: number of rows of the covariance matrix \n",
        "   * Returns: \n",
        "   *   A nrows x nrows MA1 covariance matrix \n",
        "   */ \n",
        "   matrix cov_matrix_ma1(real ma, real sigma, int nrows) { \n",
        "     matrix[nrows, nrows] mat; \n",
        "     mat <- diag_matrix(rep_vector(1 + ma^2, nrows)); \n",
        "     if (nrows > 1) { \n",
        "       mat[1, 2] <- ma; \n",
        "       for (i in 2:(nrows - 1)) { \n",
        "         mat[i, i - 1] <- ma; \n",
        "         mat[i, i + 1] <- ma; \n",
        "       } \n",
        "       mat[nrows, nrows - 1] <- ma; \n",
        "     } \n",
        "     return sigma^2 * mat; \n",
        "   } \n")
      } else {
        out$fun <- paste0(out$fun,
        "  /* compute the covariance matrix for an ARMA1 process \n",
        "   * Args: \n",
        "   *   ar: AR1 autocorrelation \n",
        "   *   ma: MA1 autocorrelation \n",
        "   *   sigma: standard deviation of the ARMA1 process \n",
        "   *   nrows: number of rows of the covariance matrix \n",
        "   * Returns: \n",
        "   *   A nrows x nrows ARMA1 covariance matrix \n",
        "   */ \n",
        "   matrix cov_matrix_arma1(real ar, real ma, real sigma, int nrows) { \n",
        "     matrix[nrows, nrows] mat; \n",
        "     vector[nrows] gamma; \n",
        "     mat <- diag_matrix(rep_vector(1 + ma^2 + 2 * ar * ma, nrows)); \n",
        "     gamma[1] <- (1 + ar * ma) * (ar + ma); \n",
        "     for (i in 2:nrows) { \n",
        "       gamma[i] <- gamma[1] * pow(ar, i - 1); \n",
        "       for (j in 1:(i - 1)) { \n",
        "         mat[i, j] <- gamma[i - j]; \n",
        "         mat[j, i] <- gamma[i - j]; \n",
        "       } \n",
        "     } \n",
        "     return sigma^2 / (1 - ar^2) * mat; \n",
        "   } \n")
      }
    } else {
      if (se) {
        stop(paste("Plese set cov = TRUE in cor_arma / cor_ar / cor_ma",
                    "when using meta-analytic standard errors"))
      }
      index <- ifelse(is_multi, "m, k", "n")
      s <- ifelse(is_multi, "  ", "")
      link_fun <- c(identity = "", log = "log", inverse = "inv")[link]
      out$transD <- paste0("  matrix[N, Karma] E;  # ARMA design matrix \n",
                           "  vector[N] e;  # residuals \n") 
      out$transC1 <- "  E <- E_pre; \n" 
      out$transC2 <- paste0(
        s,"    # calculation of ARMA effects \n",
        s,"    e[n] <- ",link_fun,"(Y[",index,"]) - eta[n]", "; \n",
        s,"    for (i in 1:Karma) { \n", 
        s,"      if (n + 1 - i > 0 && n < N && tgroup[n + 1] == tgroup[n + 1 - i]) { \n",
        s,"        E[n + 1, i] <- e[n + 1 - i]; \n",
        s,"      } \n",
        s,"    } \n")
    } 
  }
  if (Karr) {
    # autoregressive effects of the response
    out$data <- paste0(out$data,
      "  # data needed for ARR effects \n",
      "  int<lower=1> Karr; \n",
      "  matrix[N, Karr] Yarr;  # ARR design matrix \n")
    out$par <- paste0(out$par,
      "  vector[Karr] arr;  # autoregressive effects of the response \n")
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
  out <- list()
  nresp <- length(response)
  if (nresp > 1) {
    if (is.linear(family)) {
      out$data <- paste0(
        "  int<lower=1> N_trait;  # number of observations per response \n",
        "  int<lower=1> K_trait;  # number of responses \n",  
        "  int NC_trait;  # number of residual correlations \n",
        "  vector[K_trait] Y[N_trait];  # response matrix \n")
      out$par <- paste0(
        "  # parameters for multivariate linear models \n",
        "  vector<lower=0>[K_trait] sigma; \n",
        "  cholesky_factor_corr[K_trait] Lrescor; \n")
      out$loop <- c(paste0(
        "  # restructure linear predictor and add REs \n",
        "  for (m in 1:N_trait) { \n",  
        "    for (k in 1:K_trait) { \n", 
        "      int n; \n",
        "      n <- (k - 1) * N_trait + m; \n"), 
        "    } \n  } \n")
      out$prior <- paste0(
        stan_prior(class = "sigma", coef = response, prior = prior),
        stan_prior(class = "Lrescor", prior = prior))
      if (family == "gaussian") {
        out$transD <- "  cholesky_factor_cov[K_trait] LSigma; \n"
        out$transC <- paste0(
          "  # compute cholesky factor of residual covariance matrix \n",
          "  LSigma <- diag_pre_multiply(sigma, Lrescor); \n")
      } else if (family %in% c("student", "cauchy")) {
        out$transD <- "  cov_matrix[K_trait] Sigma; \n"
        out$transC <- paste0(
          "  # compute residual covariance matrix \n",
          "  Sigma <- multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor)); \n")
      }
      out$genD <- paste0(
        "  matrix[K_trait,K_trait] Rescor; \n",
        "  vector<lower=-1,upper=1>[NC_trait] rescor; \n")
      out$genC <- paste0(
        "  # take only relevant parts of residual correlation matrix \n",
        "  Rescor <- multiply_lower_tri_self_transpose(Lrescor); \n",
        collapse(ulapply(2:nresp, function(i) lapply(1:(i-1), function(j)
        paste0("  rescor[",(i-1)*(i-2)/2+j,"] <- Rescor[",j,", ",i,"]; \n")))))
    } else if (!(is.zero_inflated(family) || is.hurdle(family))) {
      stop("invalid multivariate model")
    }
  }
  out
}

stan_ordinal <- function(family, link, prior = prior_frame(), 
                         partial = FALSE, threshold = "flexible") {
  # Ordinal effects in Stan
  #
  # Args:
  #   family: the model family
  #   link: the link function
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   partial: logical; are there partial effects?
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   A vector of strings containing the ordinal effects in stan language
  out <- list()
  if (is.ordinal(family)) {
    # define Stan code similar for all ordinal models
    th <- function(k, fam = family) {
      # helper function generating stan code inside ilink(.)
      sign <- ifelse(fam %in% c("cumulative", "sratio")," - ", " + ")
      ptl <- ifelse(partial, paste0(sign, "etap[n,k]"), "") 
      if (sign == " - ") {
        out <- paste0("b_Intercept[",k,"]", ptl, " - eta[n]")
      } else {
        out <- paste0("eta[n]", ptl, " - b_Intercept[",k,"]")
      }
    }  
    ilink <- stan_ilink(link)
    type <- ifelse(family == "cumulative", "ordered", "vector")
    intercept <- paste0("  ", type, "[ncat-1] b_Intercept;  # thresholds \n")
    if (threshold == "flexible") {
      out$par <- intercept
      out$prior <- stan_prior("b_Intercept", prior = prior) 
    } else if (threshold == "equidistant") {
      out$par <- paste0("  real b_Intercept1;  # threshold 1 \n",
                        "  real", if (family == "cumulative") "<lower=0>",
                        " delta;  # distance between thresholds \n")
      out$transD <- intercept
      out$transC1 <- paste0("  # compute equidistant thresholds \n",
                            "  for (k in 1:(ncat - 1)) { \n",
                            "    b_Intercept[k] <- b_Intercept1 + (k - 1.0)*delta; \n",
                            "  } \n")
      out$prior <- paste0(stan_prior(class = "b_Intercept1", prior = prior), 
                          stan_prior(class = "delta", prior = prior))
    }
    
    # generate Stan code specific for each ordinal model
    if (!(family == "cumulative" && ilink == "inv_logit")) {
      # cumulative(logit) family already has a function in Stan itself
      out$transD <- paste0(out$transD, 
        "  vector[ncat] p[N]; \n", 
        if (family != "cumulative") 
          paste0("  vector[ncat - 1] q[N]; \n"))
      out$transC2 <- "    # compute probabilities for ordinal models \n"
      if (family == "cumulative") {
        out$transC2 <- paste0(out$transC2,
        "    p[n, 1] <- ",ilink,"(",th(1),"); \n",
        "    for (k in 2:(ncat - 1)) { \n", 
        "      p[n, k] <- ",ilink,"(",th("k"),") - ",ilink,"(",th("k - 1"),"); \n", 
        "    } \n",
        "    p[n, ncat] <- 1 - ",ilink,"(",th("ncat - 1"),"); \n")
      } else if (family %in% c("sratio", "cratio")) {
        sc <- ifelse(family == "sratio", "1 - ", "")
        out$transC2 <- paste0(out$transC2,
        "    for (k in 1:(ncat - 1)) { \n",
        "      q[n, k] <- ",sc, ilink,"(",th("k"),"); \n",
        "      p[n, k] <- 1 - q[n, k]; \n",
        "      for (kk in 1:(k - 1)) p[n, k] <- p[n, k] * q[n, kk]; \n", 
        "    } \n",
        "    p[n, ncat] <- prod(q[n]); \n")
      } else if (family == "acat") {
        if (ilink == "inv_logit") {
          out$transC2 <- paste0(out$transC2,
          "    p[n, 1] <- 1.0; \n",
          "    for (k in 1:(ncat - 1)) { \n",
          "      q[n, k] <- ",th("k"),"; \n",
          "      p[n, k + 1] <- q[n, 1]; \n",
          "      for (kk in 2:k) p[n, k + 1] <- p[n, k + 1] + q[n, kk]; \n",
          "      p[n, k + 1] <- exp(p[n, k + 1]); \n",
          "    } \n",
          "    p[n] <- p[n] / sum(p[n]); \n")
        } else {
          out$transC2 <- paste0(out$transC2,                   
          "    for (k in 1:(ncat - 1)) \n",
          "      q[n, k] <- ",ilink,"(",th("k"),"); \n",
          "    for (k in 1:ncat) { \n",     
          "      p[n, k] <- 1.0; \n",
          "      for (kk in 1:(k - 1)) p[n, k] <- p[n, k] * q[n, kk]; \n",
          "      for (kk in k:(ncat - 1)) p[n, k] <- p[n, k] * (1-q[n, kk]); \n",      
          "    } \n",
          "    p[n] <- p[n] / sum(p[n]); \n")
        }
      }
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
  out <- list()
  if (is.zero_inflated(family) || is.hurdle(family)) {
    out$data <- paste0(
      "  int<lower=1> N_trait;  # number of obs per response \n",
      "  ", ifelse(family == "hurdle_gamma", "real", "int"),
      " Y[N_trait];  # response variable \n")
    if (family == "zero_inflated_poisson") {
      out$fun <- paste0(out$fun, 
      "  /* zero-inflated poisson log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   eta: linear predictor for poisson part \n",
      "   *   eta_zi: linear predictor for zero-inflation part \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real zero_inflated_poisson_log(int y, real eta, real eta_zi) { \n",
      "     if (y == 0) { \n",
      "       return log_sum_exp(bernoulli_logit_log(1, eta_zi), \n",
      "                          bernoulli_logit_log(0, eta_zi) + \n",
      "                          poisson_log_log(0, eta)); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_zi) + \n", 
      "              poisson_log_log(y, eta); \n",
      "     } \n",
      "   } \n")
    } else if (family == "zero_inflated_negbinomial") {
      out$fun <- paste0(out$fun, 
      "  /* zero-inflated negative binomial log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   eta: linear predictor for negative binomial part \n",
      "   *   eta_zi: linear predictor for zero-inflation part \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real zero_inflated_neg_binomial_2_log(int y, real eta, real eta_zi, \n",
      "                                         real shape) { \n",
      "     if (y == 0) { \n",
      "       return log_sum_exp(bernoulli_logit_log(1, eta_zi), \n",
      "                          bernoulli_logit_log(0, eta_zi) + \n",
      "                          neg_binomial_2_log_log(0, eta, shape)); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_zi) + \n", 
      "              neg_binomial_2_log_log(y, eta, shape); \n",
      "     } \n",
      "   } \n")
    } else if (family == "zero_inflated_binomial") {
      out$fun <- paste0(out$fun, 
      "  /* zero-inflated binomial log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   eta: linear predictor for binomial part \n",
      "   *   eta_zi: linear predictor for zero-inflation part \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real zero_inflated_binomial_log(int y, int trials, real eta, \n",
      "                                         real eta_zi) { \n",
      "     if (y == 0) { \n",
      "       return log_sum_exp(bernoulli_logit_log(1, eta_zi), \n",
      "                          bernoulli_logit_log(0, eta_zi) + \n",
      "                          binomial_logit_log(0, trials, eta)); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_zi) + \n", 
      "              binomial_logit_log(y, trials, eta); \n",
      "     } \n",
      "   } \n")
    } else if (family == "hurdle_poisson") {
      out$fun <- paste0(out$fun, 
      "  /* hurdle poisson log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   eta: linear predictor for poisson part \n",
      "   *   eta_hu: linear predictor for hurdle part \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real hurdle_poisson_log(int y, real eta, real eta_hu) { \n",
      "     if (y == 0) { \n",
      "       return bernoulli_logit_log(1, eta_hu); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_hu) + \n", 
      "              poisson_log_log(y, eta) - \n",
      "              log(1 - exp(-exp(eta))); \n",
      "     } \n",
      "   } \n")
    } else if (family == "hurdle_negbinomial") {
      out$fun <- paste0(out$fun, 
      "  /* hurdle negative binomial log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   eta: linear predictor for negative binomial part \n",
      "   *   eta_hu: linear predictor for hurdle part \n",
      "   *   shape: shape parameter of negative binomial distribution \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real hurdle_neg_binomial_2_log(int y, real eta, real eta_hu, \n", 
      "                                  real shape) { \n",
      "     if (y == 0) { \n",
      "       return bernoulli_logit_log(1, eta_hu); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_hu) + \n", 
      "              neg_binomial_2_log_log(y, eta, shape) - \n",
      "              log(1 - (shape / (exp(eta) + shape))^shape); \n",
      "     } \n",
      "   } \n")
    } else if (family == "hurdle_gamma") {
      out$fun <- paste0(out$fun, 
      "  /* hurdle gamma log-PDF of a single response \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   shape: shape parameter of gamma distribution \n",
      "   *   eta: linear predictor for gamma part \n",
      "   *   eta_hu: linear predictor for hurdle part \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real hurdle_gamma_log(real y, real shape, real eta, \n", 
      "                         real eta_hu) { \n",
      "     if (y == 0) { \n",
      "       return bernoulli_logit_log(1, eta_hu); \n",
      "     } else { \n",
      "       return bernoulli_logit_log(0, eta_hu) + \n", 
      "              gamma_log(y, shape, shape / exp(eta)); \n",
      "     } \n",
      "   } \n")
    }
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
  out <- list()
  if (family == "inverse.gaussian") {
    out$data <- paste0(
      "  # quantities for the inverse gaussian distribution \n",
      "  vector[N] sqrt_Y;  # sqrt(Y) \n")
    if (weights || cens || trunc) {
      out$data <- paste0(out$data, "  vector[N] log_Y;  # log(Y) \n")
      out$fun <- paste0(out$fun,
      "  /* inverse Gaussian log-PDF for a single response (for data only) \n",
      "   * Copyright Stan Development Team 2015 \n",
      "   * Args: \n",
      "   *   y: the response value \n",
      "   *   mu: positive mean parameter \n",
      "   *   shape: positive shape parameter \n",
      "   *   log_y: precomputed log(y) \n",
      "   *   sqrt_y: precomputed sqrt(y) \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real inv_gaussian_log(real y, real mu, real shape, \n", 
      "                         real log_y, real sqrt_y) { \n",
      "     return 0.5 * log(shape / (2 * pi())) - \n", 
      "            1.5 * log_y - \n",
      "            0.5 * shape * square((y - mu) / (mu * sqrt_y)); \n",
      "   } \n")
    } else {
      out$data <- paste0(out$data, "  real log_Y;  # sum(log(Y)) \n")
      out$fun <- paste0(out$fun, 
      "  /* vectorized inverse Gaussian log-PDF (for data only) \n",
      "   * Copyright Stan Development Team 2015 \n",
      "   * Args: \n",
      "   *   y: response vector \n",
      "   *   mu: positive mean parameter vector \n",
      "   *   shape: positive shape parameter \n",
      "   *   sum_log_y: precomputed sum of log(y) \n",
      "   *   sqrt_y: precomputed sqrt(y) \n",
      "   * Returns: \n", 
      "   *   a scalar to be added to the log posterior \n",
      "   */ \n",
      "   real inv_gaussian_log(vector y, vector mu, real shape, \n", 
      "                         real sum_log_y, vector sqrt_y) { \n",
      "     return 0.5 * rows(y) * log(shape / (2 * pi())) - \n", 
      "            1.5 * sum_log_y - \n",
      "            0.5 * shape * dot_self((y - mu) ./ (mu .* sqrt_y)); \n",
      "   } \n")
    } 
    if (cens || trunc) {
      out$fun <- paste0(out$fun,
      "  /* inverse Gaussian log-CDF for a single quantile \n",
      "   * Args: \n",
      "   *   y: a quantile \n",
      "   *   mu: positive mean parameter \n",
      "   *   shape: positive shape parameter \n",
      "   *   log_y: ignored (cdf and pdf should have the same args) \n",
      "   *   sqrt_y: precomputed sqrt(y) \n",
      "   * Returns: \n",
      "   *   log(P(Y <= y)) \n",
      "   */ \n",
      "  real inv_gaussian_cdf_log(real y, real mu, real shape, \n", 
      "                            real log_y, real sqrt_y) { \n",
      "    return log(Phi(sqrt(shape) / sqrt_y * (y / mu - 1)) + \n",
      "               exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt_y * (y / mu + 1))); \n",
      "  } \n",
      "  /* inverse Gaussian log-CCDF for a single quantile \n",
      "   * Args: \n",
      "   *   y: a quantile \n",
      "   *   mu: positive mean parameter \n",
      "   *   shape: positive shape parameter \n",
      "   *   log_y: ignored (ccdf and pdf should have the same args) \n",
      "   *   sqrt_y: precomputed sqrt(y) \n",
      "   * Returns: \n",
      "   *   log(P(Y > y)) \n",
      "   */ \n",
      "  real inv_gaussian_ccdf_log(real y, real mu, real shape, \n",
      "                             real log_y, real sqrt_y) { \n",
      "    return log(1 - Phi(sqrt(shape) / sqrt_y * (y / mu - 1)) - \n",
      "               exp(2 * shape / mu) * Phi(-sqrt(shape) / sqrt_y * (y / mu + 1))); \n",
      "  } \n")
    }
  }
  out
}

stan_misc_functions <- function(link = "identity", kronecker = FALSE) {
  # stan code for user defined functions
  #
  # Args:
  #   link: the link function
  #   kronecker: logical; is the kronecker product needed?
  #
  # Returns:
  #   a string containing defined functions in stan code
  out <- NULL
  if (link == "cauchit") {
    out <- paste0(out,
    "  /* compute the inverse of the cauchit link \n",
    "   * Args: \n",
    "   *   y: the real value to be transformed \n",
    "   * Returns: \n",
    "   *   a scalar in (0,1) \n",
    "   */ \n",
    "  real inv_cauchit(real y) { \n",
    "    real p; \n",
    "    p <- cauchy_cdf(y, 0, 1); \n",
    "    return p; \n",
    "  } \n")
  }
  if (kronecker) {
    out <- paste0(out,
    "  /* calculate the cholesky factor of a kronecker covariance matrix \n",
    "   * Args: \n",
    "   *   X: a covariance matrix \n",
    "   *   L: cholesky factor of another covariance matrix \n",
    "   *   sd: standard deviations for scaling \n",
    "   * Returns: \n",
    "   *   cholesky factor of kronecker(X, L * L') \n", 
    "   */ \n",
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
    "            kron[(k-1) * rC+l, (i-1) * rC+j] <- sd[l] * sd[j] * X[k,i] * C[l,j]; \n",
    "          } \n",
    "        } \n",
    "      } \n",
    "    } \n",
    "    return cholesky_decompose(kron); \n",
    "  } \n",
    "  /* turn a vector into a 2 dimensional array \n",
    "   * Args: \n",
    "   *   X: a vector \n",
    "   *   N: first dimension of the desired array \n",
    "   *   K: second dimension of the desired array \n",
    "   * Returns: \n",
    "   *   an array of dimension N x K \n",
    "   */ \n",
    "  vector[] to_array(vector X, int N, int K) { \n",
    "    vector[K] Y[N]; \n",
    "    for (i in 1:N) \n",
    "      Y[i] <- segment(X, (i - 1) * K + 1, K); \n",
    "    return Y; \n",
    "  } \n")
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
    keep2 <- which(user_prior$group == group | !nchar(user_prior$group))
    user_prior <- user_prior[keep2, ]
  }
  if (!nchar(class) && nrow(user_prior)) {
    # increment_log_prob statements are directly put into the Stan code
    return(collapse(s, user_prior$prior, "; \n"))
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
  if (class == "b" && !is.null(attr(prior, "hs_df"))) {
    # add horseshoe shrinkage priors
    hs_shrinkage_priors <- paste0(
      "  hs_local ~ student_t(", attr(prior, "hs_df"), ", 0, 1); \n",
      "  hs_global ~ cauchy(0, 1); \n")
    out <- c(hs_shrinkage_priors, out)
  }
  return(collapse(out))
}

stan_rngprior <- function(sample.prior, prior, family = "gaussian",
                          hs_df = NULL) {
  # stan code to sample from priors seperately
  #
  # Args:
  #   sample.prior: take samples from priors?
  #   prior: the character string taken from stan_prior
  #   family: the model family
  #   hs_df: hs_df degrees of freedom
  #
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  out <- list()
  if (sample.prior) {
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
                   family == "cumulative" & grepl("^delta$", pars)
    if (any(bound)) {  
      # bounded parameters have to be sampled in the model block
      lower_bound <- ifelse(pars[bound] == "nu", 1, 0)
      out$par <- paste0("  # parameters to store prior samples \n",
                        collapse("  real<lower=", lower_bound, "> ", 
                                 "prior_", pars[bound], "; \n"))
      out$model <- paste0("  # additionally draw samples from priors \n",
                          collapse("  prior_", pars[bound] ," ~ ",
                            dis[bound], args[bound]," \n"))
    }
    if (any(!bound)) {  
      # unbounded parameters can be sampled in the generatated quantities block
      if (!is.null(hs_df)) {
        args[match("b", pars)] <- "(0, prior_hs_local * prior_hs_global);" 
      } 
      out$genD <- collapse("  real prior_", pars[!bound], "; \n")
      out$genC <- paste0("  # additionally draw samples from priors \n",
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
