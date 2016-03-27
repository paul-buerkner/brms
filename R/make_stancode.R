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
make_stancode <- function(formula, data = NULL, family = gaussian(), 
                          prior = NULL, autocor = NULL, 
                          nonlinear = NULL, partial = NULL, 
                          threshold = c("flexible", "equidistant"),
                          sparse = FALSE,  cov_ranef = NULL, 
                          sample_prior = FALSE, stan_funs = NULL, 
                          save_model = NULL, ...) {
  dots <- list(...)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef)
  sample_prior <- use_alias(sample_prior, dots$sample.prior)
  save_model <- use_alias(save_model, dots$save.model)
  dots[c("cov.ranef", "sample.prior", "save.model")] <- NULL
  # some input checks 
  if (!(is.null(data) || is.list(data)))
    stop("argument 'data' must be a data.frame or list", call. = FALSE)
  family <- check_family(family) 
  nonlinear <- nonlinear2list(nonlinear) 
  formula <- update_formula(formula, data = data, family = family, 
                            partial = partial, nonlinear = nonlinear)
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  et <- extract_time(autocor$formula)  
  ee <- extract_effects(formula, family = family, et$all, 
                        nonlinear = nonlinear)
  prior <- check_prior(prior, formula = formula, data = data, 
                       family = family, autocor = autocor, 
                       threshold = threshold, nonlinear = nonlinear) 
  data <- update_data(data, family = family, effects = ee, et$group)
  
  # flags to indicate of which type family is
  # see misc.R for details
  is_linear <- is.linear(family)
  is_ordinal <- is.ordinal(family)
  is_skewed <- is.skewed(family)
  is_count <- is.count(family)
  is_hurdle <- is.hurdle(family)
  is_zero_inflated <- is.zero_inflated(family)
  is_categorical <- is.categorical(family)
  is_multi <- is_linear && length(ee$response) > 1L
  is_forked <- is.forked(family)
  has_sigma <- has_sigma(family, autocor = autocor, se = ee$se, 
                         is_multi = is_multi)
  has_shape <- has_shape(family)
  offset <- !is.null(model.offset(data)) 
  trunc <- get_boundaries(ee$trunc)  
  add <- is.formula(ee[c("weights", "cens", "trunc")])
  
  ranef <- gather_ranef(ee, data = data, forked = is_forked)
  if (length(nonlinear)) {
    text_nonlinear <- stan_nonlinear(ee, data = data, family = family, 
                                     add = add, cov_ranef = cov_ranef,
                                     prior = prior, disp = is.formula(ee$disp))
    # currently unused objects in non-linear models
    text_fixef <- text_monef <- text_ranef <- text_csef <- 
      text_eta <- intercepts <- NULL
  } else {
    # generate fixed effects code
    intercepts <- names(get_intercepts(ee, family = family, data = data))
    fixef <- colnames(get_model_matrix(rhs(ee$fixed), data, forked = is_forked,
                                       cols2remove = intercepts))
    text_fixef <- stan_fixef(fixef = fixef, intercepts = intercepts, 
                             family = family, prior = prior, 
                             sparse = sparse, threshold = threshold)
    # generate code for category specific effects
    csef <- colnames(get_model_matrix(ee$cse, data))
    text_csef <- stan_csef(csef = csef, prior = prior)
    # generate code for monotonous effects
    monef <- colnames(get_model_matrix(ee$mono, data))
    text_monef <- stan_monef(monef, prior = prior)
    # generate random effects code
    # call stan_ranef for each random term seperately
    text_ranef <- lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
                         names_cov_ranef = names(cov_ranef), prior = prior)
    # combine random effects stan code of different grouping factors by names
    text_ranef <- collapse_lists(text_ranef)
    # generate stan code for the linear predictor
    text_eta <- stan_eta(family = family, fixef = fixef, ranef = ranef,
                         nint = length(intercepts), monef = monef, csef = csef, 
                         autocor = autocor, offset = offset, sparse = sparse,
                         add = add, disp = is.formula(ee$disp), is_multi = is_multi)
    text_nonlinear <- list()
  }
  # generate stan code for the likelihood
  text_llh <- stan_llh(family, se = is.formula(ee$se),  
                       weights = is.formula(ee$weights),
                       trials = is.formula(ee$trials),
                       cens = is.formula(ee$cens),
                       disp = is.formula(ee$disp),
                       trunc = trunc, autocor = autocor,
                       cse = is.formula(ee$cse), 
                       nresp = length(ee$response))
  # generate stan code specific to certain models
  text_arma <- stan_arma(family, autocor = autocor, prior = prior,
                         has_se = is.formula(ee$se), nonlinear = nonlinear,
                         has_disp = is.formula(ee$disp), is_multi = is_multi)
  text_multi <- stan_multi(family, response = ee$response, prior = prior)
  text_ordinal <- stan_ordinal(family, prior = prior, cse = length(csef), 
                               threshold = threshold)
  text_categorical <- stan_categorical(family)
  text_zi_hu <- stan_zero_inflated_hurdle(family)
  text_2PL <- stan_2PL(family)
  text_inv_gaussian <- stan_inv_gaussian(family, weights = is.formula(ee$weights),
                                         cens = is.formula(ee$cens),
                                         trunc = is.formula(ee$trunc))
  text_disp <- stan_disp(is.formula(ee$disp), family = family)
  kronecker <- needs_kronecker(ranef, names_cov_ranef = names(cov_ranef))
  text_misc_funs <- stan_misc_functions(family = family, kronecker = kronecker)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_fixef$prior,
    text_monef$prior,
    text_csef$prior,
    text_ordinal$prior,
    text_ranef$prior,
    text_nonlinear$prior,
    text_arma$prior,
    text_multi$prior,
    if (has_sigma) 
      stan_prior(class = "sigma", coef = ee$response, prior = prior), 
    if (has_shape) 
      stan_prior(class = "shape", prior = prior),
    if (family$family == "student") 
      stan_prior(class = "nu", prior = prior),
    if (family$family %in% c("beta", "zero_inflated_beta")) 
      stan_prior(class = "phi", prior = prior),
    stan_prior(class = "", prior = prior))
  # generate code to additionally sample from priors if sample_prior = TRUE
  text_rngprior <- stan_rngprior(sample_prior = sample_prior, 
                                 prior = text_prior, family = family,
                                 hs_df = attr(prior, "hs_df"))
  
  # generate functions block
  text_functions <- paste0(
    "// This Stan code was generated with the R package 'brms'. \n",
    "// We recommend generating the data with the 'make_standata' function. \n",
    "functions { \n",
      text_misc_funs,
      text_monef$fun,
      text_arma$fun,
      text_ordinal$fun,
      text_zi_hu$fun,
      text_inv_gaussian$fun,
      stan_funs,
    "} \n")
  
  # generate data block
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  N_bin <- ifelse(is.formula(ee$trials), "[N]", "")
  trait <- ifelse(is_multi || is_forked || is_categorical, "_trait", "")
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    if (is_multi) {
      text_multi$data
    } else if (is_categorical) {
      text_categorical$data
    } else if (is_forked) {
      paste0("  int<lower=1> N_trait;  // number of observation per trait \n",
             "  ", ifelse(use_real(family), "real", "int"),
             " Y[N_trait];  // response variable \n")
    } else if (use_real(family)) {
      "  vector[N] Y;  // response variable \n"
    } else if (use_int(family)) {
      "  int Y[N];  // response variable \n"
    },
    text_fixef$data,
    text_monef$data,
    text_csef$data,
    text_ranef$data,
    text_nonlinear$data,
    text_ordinal$data,
    text_arma$data,
    text_inv_gaussian$data,
    text_disp$data,
    if (has_trials(family))
      paste0("  int trials", N_bin, ";  // number of trials \n"),
    if (offset)
      "  vector[N] offset;  // added to the linear predictor \n",
    if (is.formula(ee$se) && !use_cov(autocor))
      "  vector<lower=0>[N] se;  // SEs for meta-analysis \n",
    if (is.formula(ee$weights))
      paste0("  vector<lower=0>[N", trait, "] weights;  // model weights \n"),
    if (is.formula(ee$cens))
      paste0("  vector[N", trait, "] cens;  // indicates censoring \n"),
    if (trunc$lb > -Inf)
      paste0("  ", ifelse(use_int(family), "int", "real"), " lb;",  
             "  // lower bound for truncation; \n"),
    if (trunc$ub < Inf)
      paste0("  ", ifelse(use_int(family), "int", "real"), " ub;",  
             "  // upper bound for truncation; \n"),
    "} \n")
  
  # generate transformed parameters block
  text_transformed_data <- paste0(
    "transformed data { \n",
       text_categorical$tdataD,
       text_fixef$tdataD,
       text_ranef$tdataD, 
       text_categorical$tdataC,
       text_fixef$tdataC,
       text_ranef$tdataC,
    "} \n")
  
  # generate parameters block
  text_parameters <- paste0(
    "parameters { \n",
    text_fixef$par,
    text_ordinal$par,
    text_monef$par,
    text_csef$par,
    text_ranef$par,
    text_nonlinear$par,
    text_arma$par,
    text_multi$par,
    if (has_sigma)
      "  real<lower=0> sigma;  // residual SD \n",
    if (family$family == "student") 
      "  real<lower=1> nu;  // degrees of freedom \n",
    if (has_shape) 
      "  real<lower=0> shape;  // shape parameter \n",
    if (family$family %in% c("beta", "zero_inflated_beta")) 
      "  real<lower=0> phi;  // precision parameter \n",
    if (!is.null(attr(prior, "hs_df"))) 
      paste0("  // horseshoe shrinkage parameters \n",
             "  vector<lower=0>[K] hs_local; \n",
             "  real<lower=0> hs_global; \n"),
    text_rngprior$par,
    "} \n")
  
  # generate transformed parameters block
  # loop over all observations in transformed parameters if necessary
  make_loop <- nrow(ee$random) || (Kar || Kma) && !use_cov(autocor) ||  
               length(intercepts) > 1L || length(ee$mono) ||
               isTRUE(text_eta$transform) || length(nonlinear)
  if (make_loop && !is_multi) {
    text_loop <- c("  for (n in 1:N) { \n", "  } \n")
  } else if (is_multi) {
    text_loop <- text_multi$loop
  } else {
    text_loop <- rep("", 2)
  }
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_eta$transD, 
      text_nonlinear$transD,
      text_disp$transD,
      text_arma$transD, 
      text_ordinal$transD,
      text_multi$transD,
      text_2PL$transD,
      text_ranef$transD, 
      text_eta$transC1,
      text_disp$transC,
      text_arma$transC1, 
      text_ordinal$transC1, 
      text_ranef$transC, 
      text_nonlinear$transC1,
      text_loop[1],
        text_eta$transC2, 
        text_nonlinear$transC2,
        text_arma$transC2, 
        text_ordinal$transC2, 
        text_eta$transC3, 
      text_loop[2],
      text_multi$transC,
      text_2PL$transC,
    "} \n")
  
  # generate model block
  lp_pre_needed <- is.formula(ee$weights) && !is.formula(ee$cens)
  text_model <- paste0(
    "model { \n",
      if (lp_pre_needed) 
        paste0("  vector[N",trait,"] lp_pre; \n"),
      "  // prior specifications \n", 
      text_prior, 
      "  // likelihood contribution \n",
      text_llh, 
      if (lp_pre_needed)  
        "  increment_log_prob(dot_product(weights, lp_pre)); \n",
      text_rngprior$model,
    "} \n")
  
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_fixef$genD,
      text_multi$genD, 
      text_ranef$genD, 
      text_nonlinear$genD,
      text_rngprior$genD,
      text_fixef$genC,
      text_multi$genC, 
      text_ranef$genC, 
      text_nonlinear$genC,
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
  
  # expand '#include' statements by calling rstan::stanc_builder
  if (!isTRUE(dots$testmode)) { 
    temp_file <- tempfile(fileext = ".stan")
    cat(complete_model, file = temp_file) 
    isystem <- system.file("chunks", package = "brms")
    complete_model <- rstan::stanc_builder(file = temp_file, isystem = isystem,
                                           obfuscate_model_name = TRUE)
    complete_model$model_name <- model_name(family)
    class(complete_model$model_code) <- c("character", "brmsmodel")
    if (is.character(save_model)) {
      cat(complete_model$model_code, file = save_model)
    }
    if (!isTRUE(dots$brm_call)) {
      complete_model <- complete_model$model_code
    }
  }
  complete_model
}