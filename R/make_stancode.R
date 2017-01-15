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
make_stancode <- function(formula, data, family = gaussian(), 
                          prior = NULL, autocor = NULL, nonlinear = NULL,
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
  family <- check_family(family)
  formula <- amend_formula(formula, data = data, family = family,
                           nonlinear = nonlinear)
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  ee <- extract_effects(formula, family = family, autocor = autocor)
  prior <- check_prior(prior, formula = formula, data = data,
                       family = family, autocor = autocor,
                       threshold = threshold, 
                       warn = !isTRUE(dots$brm_call))
  prior_only <- identical(sample_prior, "only")
  sample_prior <- if (prior_only) FALSE else sample_prior
  data <- update_data(data, family = family, effects = ee)
  
  # flags to indicate the family type
  is_categorical <- is_categorical(family)
  is_mv <- is_linear(family) && length(ee$response) > 1L
  is_forked <- is_forked(family)
  has_cens <- has_cens(ee$cens, data = data)
  bounds <- get_bounds(ee$trunc, data = data)
  
  ranef <- tidy_ranef(ee, data = data)
  if (length(ee$nlpars)) {
    text_pred <- stan_nonlinear(ee, data = data, family = family, 
                                ranef = ranef, prior = prior)
  } else {
    if (length(ee$response) > 1L) {
      text_pred <- stan_effects_mv(ee, data = data, family = family,
                                   prior = prior, ranef = ranef, 
                                   autocor = autocor, sparse = sparse)
    } else {
      text_pred <- stan_effects(ee, data = data, family = family,
                                prior = prior, ranef = ranef, 
                                autocor = autocor, sparse = sparse,
                                threshold = threshold)
    }
  }
  text_auxpars <- stan_auxpars(ee, data = data, family = family, 
                               ranef = ranef, prior = prior, 
                               autocor = autocor)
  # due to the ID syntax, group-level effects are evaluated separately
  text_ranef <- lapply(unique(ranef$id), stan_ranef, ranef = ranef, 
                       prior = prior, cov_ranef = cov_ranef)
  text_ranef <- collapse_lists(text_ranef)
  
  # generate Stan code of the likelihood
  text_llh <- stan_llh(family, effects = ee, data = data, autocor = autocor)
  # generate Stan code specific to certain models
  text_autocor <- stan_autocor(autocor, effects = ee, family = family,
                               prior = prior)
  text_mv <- stan_mv(family, response = ee$response, prior = prior)
  disc <- "disc" %in% names(ee$auxpars) || isTRUE(ee$fauxpars$disc != 1)
  text_ordinal <- stan_ordinal(family, prior = prior, cs = has_cs(ee), 
                               disc = disc, threshold = threshold)
  text_families <- stan_families(family)
  text_se <- stan_se(is.formula(ee$se))
  text_cens <- stan_cens(has_cens, family = family)
  text_disp <- stan_disp(ee, family = family)
  kronecker <- stan_needs_kronecker(ranef, names_cov_ranef = names(cov_ranef))
  text_misc_funs <- stan_misc_functions(family = family, kronecker = kronecker)
  text_monotonic <- stan_monotonic(text_pred)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_pred$prior,
    text_auxpars$prior,
    text_ranef$prior,
    text_ordinal$prior,
    text_autocor$prior,
    text_mv$prior,
    stan_prior(class = "", prior = prior))
  
  # generate functions block
  text_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions { \n",
      text_misc_funs,
      text_monotonic,
      text_autocor$fun,
      text_ordinal$fun,
      text_families$fun,
      stan_funs,
    "} \n")
  
  # generate data block
  rtype <- ifelse(use_int(family), "int", "real")
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    if (is_mv) {
      text_mv$data
    } else if (use_real(family)) {
      # don't use real Y[n]
      "  vector[N] Y;  // response variable \n"
    } else if (use_int(family)) {
      "  int Y[N];  // response variable \n"
    },
    text_pred$data,
    text_auxpars$data,
    text_ranef$data,
    text_ordinal$data,
    text_families$data,
    text_autocor$data,
    text_cens$data,
    text_disp$data,
    text_se$data,
    if (has_trials(family))
      "  int trials[N];  // number of trials \n",
    if (is.formula(ee$weights))
      "  vector<lower=0>[N] weights;  // model weights \n",
    if (is.formula(ee$dec))
      "  int<lower=0,upper=1> dec[N];  // decisions \n",
    if (any(bounds$lb > -Inf))
      paste0("  ", rtype, " lb[N];  // lower truncation bounds; \n"),
    if (any(bounds$ub < Inf))
      paste0("  ", rtype, " ub[N];  // upper truncation bounds \n"),
    "  int prior_only;  // should the likelihood be ignored? \n",
    "} \n")
  
  # generate transformed parameters block
  text_transformed_data <- paste0(
    "transformed data { \n",
       text_families$tdataD,
       text_pred$tdataD,
       text_auxpars$tdataD,
       text_se$tdataD,
       text_autocor$tdataD,
       text_families$tdataC,
       text_pred$tdataC,
       text_auxpars$tdataC,
       text_se$tdataC,
       text_autocor$tdataC,
    "} \n")
  
  # generate parameters block
  text_parameters <- paste0(
    text_pred$par,
    text_auxpars$par,
    text_ranef$par,
    text_ordinal$par,
    text_autocor$par,
    text_mv$par)
  text_rngprior <- stan_rngprior(sample_prior = sample_prior, 
                                 par_declars = text_parameters,
                                 prior = text_prior, family = family,
                                 prior_attr = attributes(prior))
  text_parameters <- paste0(
    "parameters { \n",
      text_parameters,
      text_rngprior$par,
    "} \n")
  
  # generate transformed parameters block
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_pred$transD,
      text_auxpars$transD,
      text_ranef$transD,
      text_autocor$transD, 
      text_ordinal$transD,
      text_mv$transD,
      text_pred$transC1,
      text_auxpars$transC1,
      text_ranef$transC1,
      text_autocor$transC1, 
      text_ordinal$transC1, 
      text_mv$transC1,
    "} \n")
  
  # generate model block
  # list auxpars before pred as part of fixing issue #124
  text_model_loop <- paste0(
    text_auxpars$modelC2, 
    text_pred$modelC2, 
    text_autocor$modelC2,
    text_auxpars$modelC3,
    text_pred$modelC3
  )
  if (isTRUE(nzchar(text_model_loop))) {
    text_model_loop <- paste0(
      "  for (n in 1:N) { \n", text_model_loop, "  } \n"
    )
  }
  text_lp_pre <- list()
  if (is.formula(ee$weights) && !is.formula(ee$cens)) {
    text_lp_pre <- list(
      modelD = "  vector[N] lp_pre; \n",
      modelC = "    target += dot_product(weights, lp_pre); \n"
    )
  }
  text_model <- paste0(
    "model { \n",
      text_pred$modelD,
      text_auxpars$modelD,
      text_disp$modelD,
      text_autocor$modelD,
      text_lp_pre$modelD,
      text_auxpars$modelC1,
      text_pred$modelC1,
      text_autocor$modelC1, 
      text_disp$modelC1,
      text_model_loop,
      "  // prior specifications \n", 
      text_prior, 
      "  // likelihood contribution \n",
      "  if (!prior_only) { \n  ",
      text_llh, 
      text_lp_pre$modelC,
      "  } \n", 
      text_rngprior$model,
    "} \n")
  
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_pred$genD,
      text_auxpars$genD,
      text_ranef$genD,
      text_mv$genD, 
      text_rngprior$genD,
      text_pred$genC,
      text_auxpars$genC,
      text_ranef$genC,
      text_mv$genC, 
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