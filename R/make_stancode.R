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
                          threshold = c("flexible", "equidistant"),
                          sparse = FALSE,  cov_ranef = NULL, 
                          sample_prior = c("no", "yes", "only"), 
                          stan_funs = NULL, save_model = NULL, 
                          silent = FALSE, ...) {
  dots <- list(...)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots[["cov.ranef"]])
  sample_prior <- use_alias(sample_prior, dots[["sample.prior"]])
  save_model <- use_alias(save_model, dots[["save.model"]])
  dots[c("cov.ranef", "sample.prior", "save.model")] <- NULL
  # some input checks
  formula <- validate_formula(
    formula, data = data, family = family, 
    autocor = autocor, threshold = threshold
  )
  bterms <- parse_bf(formula)
  sample_prior <- check_sample_prior(sample_prior)
  prior <- check_prior(
    prior, formula, data = data, sample_prior = sample_prior,
    warn = !isTRUE(dots$brm_call)
  )
  data <- update_data(data, bterms = bterms)
  ranef <- tidy_ranef(bterms, data = data)
  
  text_effects <- stan_effects(
    bterms, data = data, ranef = ranef, 
    prior = prior, sparse = sparse
  )
  # the ID syntax requires group-level effects to be evaluated separately
  text_ranef <- collapse_lists(ls = lapply(
    X = unique(ranef$id), FUN = stan_re,
    ranef = ranef, prior = prior, cov_ranef = cov_ranef
  ))
  text_llh <- stan_llh(bterms, data = data)
  text_global_defs <- stan_global_defs(
    bterms, prior = prior, ranef = ranef, cov_ranef = cov_ranef
  )
  text_Xme <- stan_Xme(bterms, prior = prior)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_effects$prior,
    text_ranef$prior,
    text_Xme$prior,
    stan_prior(class = "", prior = prior)
  )
  # generate functions block
  text_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions { \n",
      text_global_defs$fun,
      stan_funs,
    "} \n"
  )
  # generate data block
  text_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    text_effects$data,
    text_ranef$data,
    text_Xme$data,
    "  int prior_only;  // should the likelihood be ignored? \n",
    "} \n"
  )
  # generate transformed parameters block
  text_transformed_data <- paste0(
    "transformed data { \n",
       text_global_defs$tdataD,
       text_effects$tdataD,
       text_effects$tdataC,
    "} \n"
  )
  # generate parameters block
  text_parameters <- paste0(
    text_effects$par,
    text_ranef$par,
    text_Xme$par
  )
  text_rngprior <- stan_rngprior(
    sample_prior = sample_prior, 
    par_declars = text_parameters,
    gen_quantities = text_effects$genD,
    prior = text_prior,
    prior_special = attr(prior, "special")
  )
  text_parameters <- paste0(
    "parameters { \n",
      text_parameters,
      text_rngprior$par,
    "} \n"
  )
  # generate transformed parameters block
  text_transformed_parameters <- paste0(
    "transformed parameters { \n",
      text_effects$tparD,
      text_ranef$tparD,
      text_effects$tparC1,
      text_ranef$tparC1,
    "} \n"
  )
  # generate model block
  text_model_loop <- paste0(
    text_effects$modelC2, 
    text_effects$modelC3,
    text_effects$modelC4
  )
  if (isTRUE(nzchar(text_model_loop))) {
    text_model_loop <- paste0(
      "  for (n in 1:N) { \n", text_model_loop, "  } \n"
    )
  }
  text_model <- paste0(
    "model { \n",
      text_effects$modelD,
      text_effects$modelC1,
      text_effects$modelCgp1,
      text_model_loop,
      text_effects$modelC5,
      "  // priors including all constants \n", 
      text_prior, 
      "  // likelihood including all constants \n",
      "  if (!prior_only) { \n",
      text_llh, 
      "  } \n", 
      text_rngprior$model,
    "} \n"
  )
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_effects$genD,
      text_ranef$genD,
      text_rngprior$genD,
      text_effects$genC,
      text_ranef$genC,
      text_rngprior$genC,
    "} \n"
  )
  # combine all elements into a complete Stan model
  complete_model <- paste0(
    text_functions,
    text_data, 
    text_transformed_data, 
    text_parameters,
    text_transformed_parameters,
    text_model,
    text_generated_quantities
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
    complete_model$model_name <- name_model(family)
    class(complete_model$model_code) <- c("character", "brmsmodel")
    if (is.character(save_model)) {
      cat(complete_model$model_code, file = save_model)
    }
    if (!isTRUE(dots$brm_call)) {
      complete_model <- complete_model$model_code
    }
  } else {
    class(complete_model) <- c("character", "brmsmodel")
  }
  complete_model
}
