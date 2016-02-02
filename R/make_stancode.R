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
                          cov_ranef = NULL, sample_prior = FALSE, 
                          save_model = NULL, ...) {
  dots <- list(...)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef)
  sample_prior <- use_alias(sample_prior, dots$sample.prior)
  save_model <- use_alias(save_model, dots$save.model)
  dots[c("cov.ranef", "sample.prior", "save.model")] <- NULL
  # some input checks 
  formula <- update_formula(formula, data = data)
  family <- check_family(family) 
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  prior <- check_prior(prior, formula = formula, data = data, 
                       family = family, autocor = autocor, 
                       partial = partial, threshold = threshold,
                       nonlinear = nonlinear) 
  et <- extract_time(autocor$formula)  
  ee <- extract_effects(formula, family = family, partial, et$all, 
                        nonlinear = nonlinear)
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
  is_multi <- is_linear && length(ee$response) > 1
  is_forked <- is.forked(family)
  has_sigma <- has_sigma(family, autocor = autocor, se = ee$se, 
                         is_multi = is_multi)
  has_shape <- has_shape(family)
  offset <- !is.null(model.offset(data)) 
  trunc <- get_boundaries(ee$trunc)  
  add <- is.formula(ee[c("weights", "cens", "trunc")])
  
  ranef <- gather_ranef(ee, data = data, is_forked = is_forked)
  if (length(nonlinear)) {
    text_nonlinear <- stan_nonlinear(ee, data = data, family = family, 
                                     add = add, cov_ranef = cov_ranef,
                                     prior = prior)
    text_fixef <- text_ranef <- text_eta <- list()
  } else {
    # generate fixed effects code
    rm_intercept <- isTRUE(attr(ee$fixed, "rsv_intercept"))
    if (is_categorical) {
      X <- data.frame()
      fixef <- colnames(X)
      Xp <- get_model_matrix(ee$fixed, data, rm_intercept = rm_intercept)
      temp_list <- check_intercept(colnames(Xp))
      paref <- temp_list$names
    } else {
      X <- get_model_matrix(ee$fixed, data, is_forked = is_forked,
                            rm_intercept = rm_intercept)
      temp_list <- check_intercept(colnames(X))
      fixef <- temp_list$names
      Xp <- get_model_matrix(partial, data, rm_intercept = TRUE)
      paref <- colnames(Xp)
    }
    has_intercept <- temp_list$has_intercept
    text_fixef <- stan_fixef(fixef = fixef, paref = paref, family = family, 
                             prior = prior, threshold = threshold,
                             has_intercept = has_intercept)
    
    # generate random effects code
    # call stan_ranef for each random term seperately
    text_ranef <- lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
                         names_cov_ranef = names(cov_ranef), prior = prior)
    # combine random effects stan code of different grouping factors by names
    text_ranef <- collapse_lists(text_ranef)
    # generate stan code for the linear predictor
    text_eta <- stan_eta(family = family, fixef = fixef, ranef = ranef,
                         has_intercept = has_intercept, paref = paref, 
                         autocor = autocor, offset = offset, 
                         is_multi = is_multi, add = add)
    text_nonlinear <- list()
  }
  # generate stan code for the likelihood
  text_llh <- stan_llh(family, se = is.formula(ee$se),  
                       weights = is.formula(ee$weights),
                       trials = is.formula(ee$trials),
                       cens = is.formula(ee$cens), 
                       trunc = trunc, autocor = autocor,
                       partial = is.formula(partial),
                       is_multi = is_multi)
  trait <- ifelse(is_multi || is.forked(family), "_trait", "")
  if (is.formula(ee$cens) || is.formula(ee$weights) || is.formula(ee$trunc) ||
      is_ordinal || is_categorical || is_hurdle || is_zero_inflated) {
    text_llh <- paste0("  for (n in 1:N",trait,") { \n  ",text_llh,"  } \n")
  }
  
  # generate stan code specific to certain models
  text_arma <- stan_arma(family = family, autocor = autocor, prior = prior,
                         se = is.formula(ee$se), is_multi = is_multi)
  text_multi <- stan_multi(family = family, response = ee$response,
                           prior = prior)
  text_ordinal <- stan_ordinal(family = family, prior = prior, 
                               partial = length(paref), 
                               threshold = threshold)  
  text_zi_hu <- stan_zero_inflated_hurdle(family)
  text_2PL <- stan_2PL(family)
  text_inv_gaussian <- stan_inv_gaussian(family = family, 
                                         weights = is.formula(ee$weights),
                                         cens = is.formula(ee$cens),
                                         trunc = is.formula(ee$trunc))
  kronecker <- needs_kronecker(ranef, names_cov_ranef = names(cov_ranef))
  text_misc_funs <- stan_misc_functions(family = family, kronecker = kronecker)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_fixef$prior,
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
    "functions { \n",
      text_misc_funs,
      text_arma$fun,
      text_ordinal$fun,
      text_zi_hu$fun,
      text_inv_gaussian$fun,
    "} \n")
  
  # generate data block
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  N_bin <- ifelse(is.formula(ee$trials), "[N]", "")
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // number of observations \n", 
    if (is_multi) {
      text_multi$data
    } else if (is_forked) {
      paste0("  int<lower=1> N_trait;  // number of obs / 2 \n",
             "  ", ifelse(use_real(family), "real", "int"),
             " Y[N_trait];  // response variable \n")
    } else if (use_real(family)) {
      "  vector[N] Y;  // response variable \n"
    } else if (use_int(family)) {
      "  int Y[N];  // response variable \n"
    },
    text_fixef$data,
    text_ranef$data,
    text_nonlinear$data,
    text_arma$data,
    text_inv_gaussian$data,
    if (has_trials(family))
      paste0("  int trials", N_bin, ";  // number of trials \n"),
    if (has_cat(family))
      paste0("  int ncat;  // number of categories \n"),
    if (offset)
      "  vector[N] offset;  // added to the linear predictor \n",
    if (is.formula(ee$se) && !(use_cov(autocor) && (Kar || Kma)))
      "  vector<lower=0>[N] se;  // SEs for meta-analysis \n",
    if (is.formula(ee$weights))
      paste0("  vector<lower=0>[N",trait,"] weights;  // model weights \n"),
    if (is.formula(ee$cens))
      paste0("  vector[N",trait,"] cens;  // indicates censoring \n"),
    if (trunc$lb > -Inf)
      paste0("  ", ifelse(use_int(family), "int", "real"), " lb;",  
             "  // lower bound for truncation; \n"),
    if (trunc$ub < Inf)
      paste0("  ", ifelse(use_int(family), "int", "real"), " ub;",  
             "  // upper bound for truncation; \n"),
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
               isTRUE(text_eta$transform) || length(nonlinear)
  if (make_loop && !is_multi) {
    text_loop <- c(paste0(
      "  // if available add REs to linear predictor \n",
      "  for (n in 1:N) { \n"), "  } \n")
  } else if (is_multi) {
    text_loop <- text_multi$loop
  } else {
    text_loop <- rep("", 2)
  }
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_eta$transD, 
      text_nonlinear$transD,
      text_arma$transD, 
      text_ordinal$transD,
      text_multi$transD,
      text_2PL$transD,
      text_ranef$transD, 
      text_eta$transC1, 
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
  
  # write the stan code to a file if save_model is a character string
  class(complete_model) <- c("character", "brmsmodel")
  if (is.character(save_model)) {
    sink(save_model)
    cat(complete_model)
    sink()
  }
  complete_model
}