#' Stan Code for \pkg{brms} Models
#' 
#' Generate Stan code for \pkg{brms} models
#' 
#' @inheritParams brm
#' @param silent logical; If \code{TRUE}, warnings of
#'   the Stan parser will be suppressed.
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
                          prior = NULL, autocor = NULL,
                          cov_ranef = NULL, sparse = FALSE, 
                          sample_prior = c("no", "yes", "only"), 
                          stan_vars = NULL, stan_funs = NULL, 
                          save_model = NULL, silent = FALSE, ...) {
  dots <- list(...)
  # some input checks
  if (is.brmsfit(formula)) {
    stop2("Use 'stancode' to extract Stan code from 'brmsfit' objects.")
  }
  if (length(stan_funs) > 0) {
    stan_funs <- as_one_character(stan_funs) 
  }
  formula <- validate_formula(
    formula, data = data, family = family, autocor = autocor
  )
  bterms <- parse_bf(formula)
  sample_prior <- check_sample_prior(sample_prior)
  prior <- check_prior(
    prior, formula, data = data, sample_prior = sample_prior,
    warn = !isTRUE(dots$brm_call)
  )
  data <- update_data(data, bterms = bterms)
  ranef <- tidy_ranef(bterms, data = data)
  meef <- tidy_meef(bterms, data = data)
  stan_vars <- validate_stanvars(stan_vars)
  
  scode_effects <- stan_effects(
    bterms, data = data, prior = prior, 
    ranef = ranef, meef = meef, sparse = sparse
  )
  # the ID syntax requires group-level effects to be evaluated separately
  scode_ranef <- collapse_lists(ls = lapply(
    X = unique(ranef$id), FUN = stan_re,
    ranef = ranef, prior = prior, cov_ranef = cov_ranef
  ))
  scode_llh <- stan_llh(bterms, data = data)
  scode_global_defs <- stan_global_defs(
    bterms, prior = prior, ranef = ranef, cov_ranef = cov_ranef
  )
  scode_Xme <- stan_Xme(meef, prior = prior)
    
  # get priors for all parameters in the model
  scode_prior <- paste0(
    scode_effects$prior,
    scode_ranef$prior,
    scode_Xme$prior,
    stan_prior(class = "", prior = prior)
  )
  # generate functions block
  scode_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions { \n",
      scode_global_defs$fun,
      stan_funs,
    "} \n"
  )
  # generate data block
  scode_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    scode_effects$data,
    scode_ranef$data,
    scode_Xme$data,
    "  int prior_only;  // should the likelihood be ignored? \n",
    collapse_stanvars(stan_vars),
    "} \n"
  )
  # generate transformed parameters block
  scode_transformed_data <- paste0(
    "transformed data { \n",
       scode_global_defs$tdataD,
       scode_effects$tdataD,
       scode_effects$tdataC,
    "} \n"
  )
  # generate parameters block
  scode_parameters <- paste0(
    scode_effects$par,
    scode_ranef$par,
    scode_Xme$par
  )
  scode_rngprior <- stan_rngprior(
    sample_prior = sample_prior, 
    par_declars = scode_parameters,
    gen_quantities = scode_effects$genD,
    prior = scode_prior,
    prior_special = attr(prior, "special")
  )
  scode_parameters <- paste0(
    "parameters { \n",
      scode_parameters,
      scode_rngprior$par,
    "} \n"
  )
  # generate transformed parameters block
  scode_transformed_parameters <- paste0(
    "transformed parameters { \n",
      scode_effects$tparD,
      scode_ranef$tparD,
      scode_Xme$tparD,
      scode_effects$tparC1,
      scode_ranef$tparC1,
    "} \n"
  )
  # generate model block
  scode_model_loop <- paste0(
    scode_effects$modelC2, 
    scode_effects$modelC3,
    scode_effects$modelC4
  )
  if (isTRUE(nzchar(scode_model_loop))) {
    scode_model_loop <- paste0(
      "  for (n in 1:N) { \n", scode_model_loop, "  } \n"
    )
  }
  scode_model <- paste0(
    "model { \n",
      scode_effects$modelD,
      scode_effects$modelC1,
      scode_effects$modelCgp1,
      scode_model_loop,
      scode_effects$modelC5,
      "  // priors including all constants \n", 
      scode_prior, 
      "  // likelihood including all constants \n",
      "  if (!prior_only) { \n",
      scode_llh, 
      "  } \n", 
      scode_rngprior$model,
    "} \n"
  )
  # generate generated quantities block
  scode_generated_quantities <- paste0(
    "generated quantities { \n",
      scode_effects$genD,
      scode_ranef$genD,
      scode_Xme$genD,
      scode_rngprior$genD,
      scode_effects$genC,
      scode_ranef$genC,
      scode_rngprior$genC,
      scode_Xme$genC,
    "} \n"
  )
  # combine all elements into a complete Stan model
  complete_model <- paste0(
    scode_functions,
    scode_data, 
    scode_transformed_data, 
    scode_parameters,
    scode_transformed_parameters,
    scode_model,
    scode_generated_quantities
  )
  
  # expand '#include' statements by calling rstan::stanc_builder
  if (!isTRUE(dots$testmode)) { 
    temp_file <- tempfile(fileext = ".stan")
    cat(complete_model, file = temp_file) 
    isystem <- system.file("chunks", package = "brms")
    complete_model <- eval_silent(
      rstan::stanc_builder(
        file = temp_file, isystem = isystem, 
        obfuscate_model_name = TRUE
      ),
      type = "message", silent = silent
    )
    model_name <- paste(summarise_families(formula), "brms-model")
    complete_model$model_name <- model_name
    if (is.character(save_model)) {
      str_add(complete_model$model_code) <- "\n"
      cat(complete_model$model_code, file = save_model)
    }
    class(complete_model$model_code) <- c("character", "brmsmodel")
    if (!isTRUE(dots$brm_call)) {
      complete_model <- complete_model$model_code
    }
  } else {
    class(complete_model) <- c("character", "brmsmodel")
  }
  complete_model
}
