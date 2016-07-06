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
  formula <- update_formula(formula, data = data, family = family,
                            partial = partial, nonlinear = nonlinear)
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  et <- extract_time(autocor$formula)  
  ee <- extract_effects(formula, family = family, et$all)
  prior <- check_prior(prior, formula = formula, data = data, 
                       family = family, autocor = autocor, 
                       threshold = threshold)
  prior_only <- identical(sample_prior, "only")
  sample_prior <- if (prior_only) FALSE else sample_prior
  data <- update_data(data, family = family, effects = ee, et$group)
  
  # flags to indicate the family type
  is_categorical <- is.categorical(family)
  is_multi <- is.linear(family) && length(ee$response) > 1L
  is_forked <- is.forked(family)
  trunc <- get_boundaries(ee$trunc)
  
  intercepts <- names(get_intercepts(ee, family = family, data = data))
  if (length(ee$nonlinear)) {
    text_pred <- stan_nonlinear(ee, data = data, family = family, 
                                prior = prior, cov_ranef = cov_ranef)
  } else {
    text_pred <- stan_linear(ee, data = data, family = family, 
                             prior = prior, intercepts = intercepts, 
                             autocor = autocor, threshold = threshold, 
                             sparse = sparse, cov_ranef = cov_ranef)
  }
  text_auxpars <- stan_auxpars(ee, data = data, family = family, 
                               prior = prior, autocor = autocor, 
                               cov_ranef = cov_ranef)
  text_pred <- collapse_lists(list(text_pred, text_auxpars))
  
  # generate Stan code of the likelihood
  text_llh <- stan_llh(family, effects = ee, autocor = autocor)
  # generate Stan code specific to certain models
  text_autocor <- stan_autocor(autocor, effects = ee, family = family,
                               prior = prior)
  text_multi <- stan_multi(family, response = ee$response, prior = prior)
  text_ordinal <- stan_ordinal(family, prior = prior, cse = is.formula(ee$cse), 
                               threshold = threshold)
  text_categorical <- stan_categorical(family)
  text_forked <- stan_forked(family)
  text_inv_gaussian <- stan_inv_gaussian(family, weights = is.formula(ee$weights),
                                         cens = is.formula(ee$cens),
                                         trunc = is.formula(ee$trunc))
  text_disp <- stan_disp(ee, family = family)
  ranef <- gather_ranef(ee, data = data)
  kronecker <- stan_needs_kronecker(ranef, names_cov_ranef = names(cov_ranef))
  text_misc_funs <- stan_misc_functions(family = family, kronecker = kronecker)
  text_monotonic <- stan_monotonic(text_pred)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_pred$prior,
    text_ordinal$prior,
    text_autocor$prior,
    text_multi$prior,
    stan_prior(class = "", prior = prior))
  
  # generate functions block
  text_functions <- paste0(
    "// This Stan code was generated with the R package 'brms'. \n",
    "// We recommend generating the data with the 'make_standata' function. \n",
    "functions { \n",
      text_misc_funs,
      text_monotonic,
      text_autocor$fun,
      text_ordinal$fun,
      text_forked$fun,
      text_inv_gaussian$fun,
      stan_funs,
    "} \n")
  
  # generate data block
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  Nbin <- ifelse(is.formula(ee$trials), "[N]", "")
  trait <- ifelse(is_multi || is_forked || is_categorical, "_trait", "")
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    if (is_multi) {
      text_multi$data
    } else if (is_categorical) {
      text_categorical$data
    } else if (is_forked) {
      text_forked$data
    } else if (use_real(family)) {
      "  vector[N] Y;  // response variable \n"
    } else if (use_int(family)) {
      "  int Y[N];  // response variable \n"
    },
    text_pred$data,
    text_ordinal$data,
    text_autocor$data,
    text_inv_gaussian$data,
    text_disp$data,
    if (has_trials(family))
      paste0("  int trials", Nbin, ";  // number of trials \n"),
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
    "  int prior_only;  // should the likelihood be ignored? \n",
    "} \n")
  
  # generate transformed parameters block
  text_transformed_data <- paste0(
    "transformed data { \n",
       text_categorical$tdataD,
       text_pred$tdataD,
       text_autocor$tdataD,
       text_categorical$tdataC,
       text_pred$tdataC,
       text_autocor$tdataC,
    "} \n")
  
  # generate parameters block
  text_parameters <- paste0(
    text_pred$par,
    text_ordinal$par,
    text_autocor$par,
    text_multi$par)
  text_rngprior <- stan_rngprior(sample_prior = sample_prior, 
                                 par_declars = text_parameters,
                                 prior = text_prior, family = family,
                                 hs_df = attr(prior, "hs_df"))
  text_parameters <- paste0(
    "parameters { \n",
      text_parameters,
      text_rngprior$par,
    "} \n")
  
  # generate transformed parameters block
  # loop over all observations in transformed parameters if necessary
  make_loop <- any(nzchar(c(text_pred$transC2, text_autocor$transC2, 
                            text_ordinal$transC2, text_pred$transC3)))
  if (make_loop && !is_multi) {
    text_loop <- c("  for (n in 1:N) { \n", "  } \n")
  } else if (is_multi) {
    text_loop <- text_multi$loop
  } else {
    text_loop <- rep("", 2)
  }
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_pred$transD,
      text_disp$transD,
      text_autocor$transD, 
      text_ordinal$transD,
      text_multi$transD,
      text_forked$transD,
      text_pred$transC1,
      text_disp$transC1,
      text_autocor$transC1, 
      text_ordinal$transC1, 
      text_loop[1],
        text_pred$transC2, 
        text_autocor$transC2, 
        text_ordinal$transC2, 
        text_pred$transC3,
      text_loop[2],
      text_pred$transC4,
      text_multi$transC1,
      text_forked$transC1,
    "} \n")
  
  # generate model block
  needs_lp_pre <- is.formula(ee$weights) && !is.formula(ee$cens)
  text_model <- paste0(
    "model { \n",
      if (needs_lp_pre) 
        paste0("  vector[N", trait,"] lp_pre; \n"),
      "  // prior specifications \n", 
      text_prior, 
      "  // likelihood contribution \n",
      "  if (!prior_only) { \n  ",
      text_llh, 
      if (needs_lp_pre)
        "    target += dot_product(weights, lp_pre); \n",
      "  } \n", 
      text_rngprior$model,
    "} \n")
  
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_pred$genD,
      text_multi$genD, 
      text_rngprior$genD,
      text_pred$genC,
      text_multi$genC, 
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