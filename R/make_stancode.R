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
                          stanvars = NULL, stan_funs = NULL, 
                          save_model = NULL, silent = FALSE, ...) {
  dots <- list(...)
  silent <- as_one_logical(silent)
  # some input checks
  if (is.brmsfit(formula)) {
    stop2("Use 'stancode' to extract Stan code from 'brmsfit' objects.")
  }
  if (length(stan_funs) > 0) {
    warning2("Argument 'stan_funs' is deprecated. Please use argument ", 
             "'stanvars' instead. See ?stanvar for more help.")
    stan_funs <- as_one_character(stan_funs) 
  }
  formula <- validate_formula(
    formula, data = data, family = family, autocor = autocor
  )
  bterms <- parse_bf(formula)
  sample_prior <- check_sample_prior(sample_prior)
  prior <- check_prior(
    prior, formula, data = data, sparse = sparse, 
    sample_prior = sample_prior, warn = TRUE
  )
  data <- update_data(data, bterms = bterms)
  ranef <- tidy_ranef(bterms, data = data)
  meef <- tidy_meef(bterms, data = data)
  stanvars <- validate_stanvars(stanvars)
  
  scode_predictor <- stan_predictor(
    bterms, data = data, prior = prior, 
    ranef = ranef, meef = meef, sparse = sparse,
    stanvars = stanvars
  )
  scode_ranef <- stan_re(ranef, prior = prior, cov_ranef = cov_ranef)
  scode_llh <- stan_llh(bterms, data = data)
  scode_global_defs <- stan_global_defs(
    bterms, prior = prior, ranef = ranef, cov_ranef = cov_ranef
  )
  scode_Xme <- stan_Xme(meef, prior = prior)
    
  # get priors for all parameters in the model
  scode_prior <- paste0(
    scode_predictor$prior,
    scode_ranef$prior,
    scode_Xme$prior,
    stan_prior(class = "", prior = prior)
  )
  # generate functions block
  scode_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions { \n",
      scode_global_defs$fun,
      collapse_stanvars(stanvars, block = "functions"),
      stan_funs,
    "} \n"
  )
  # generate data block
  scode_data <- paste0(
    "data { \n",
    "  int<lower=1> N;  // total number of observations \n", 
    scode_predictor$data,
    scode_ranef$data,
    scode_Xme$data,
    "  int prior_only;  // should the likelihood be ignored? \n",
    collapse_stanvars(stanvars, block = "data"),
    "} \n"
  )
  # generate transformed parameters block
  scode_transformed_data <- paste0(
    "transformed data { \n",
       scode_global_defs$tdataD,
       scode_predictor$tdataD,
       collapse_stanvars(stanvars, block = "tdata"),
       scode_predictor$tdataC,
    "} \n"
  )
  # generate parameters block
  scode_parameters <- paste0(
    scode_predictor$par,
    scode_ranef$par,
    scode_Xme$par
  )
  scode_rngprior <- stan_rngprior(
    sample_prior = sample_prior, 
    par_declars = scode_parameters,
    gen_quantities = scode_predictor$genD,
    prior = scode_prior,
    prior_special = attr(prior, "special")
  )
  scode_parameters <- paste0(
    "parameters { \n",
      scode_parameters,
      scode_rngprior$par,
      collapse_stanvars(stanvars, block = "parameters"),
    "} \n"
  )
  # generate transformed parameters block
  scode_transformed_parameters <- paste0(
    "transformed parameters { \n",
      scode_predictor$tparD,
      scode_ranef$tparD,
      scode_Xme$tparD,
      collapse_stanvars(stanvars, block = "tparameters"),
      scode_predictor$tparC1,
      scode_ranef$tparC1,
    "} \n"
  )
  # generate model block
  scode_model_loop <- paste0(
    scode_predictor$modelC2, 
    scode_predictor$modelC3,
    scode_predictor$modelC4
  )
  if (isTRUE(nzchar(scode_model_loop))) {
    scode_model_loop <- paste0(
      "  for (n in 1:N) { \n", scode_model_loop, "  } \n"
    )
  }
  scode_model <- paste0(
    "model { \n",
      scode_predictor$modelD,
      collapse_stanvars(stanvars, block = "model"),
      scode_predictor$modelC1,
      scode_predictor$modelCgp1,
      scode_model_loop,
      scode_predictor$modelC5,
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
      scode_predictor$genD,
      scode_ranef$genD,
      scode_Xme$genD,
      scode_rngprior$genD,
      collapse_stanvars(stanvars, block = "genquant"),
      scode_predictor$genC,
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
  
  if (!isTRUE(dots$testmode)) { 
    # expand '#include' statements by calling rstan::stanc_builder
    temp_file <- tempfile(fileext = ".stan")
    cat(complete_model, file = temp_file) 
    isystem <- system.file("chunks", package = "brms")
    # get rid of diagnostic messages from parser
    complete_model <- eval_silent(
      rstan::stanc_builder(
        file = temp_file, isystem = isystem,
        obfuscate_model_name = TRUE
      ),
      type = "message", try = TRUE
    )
    complete_model <- complete_model$model_code
    str_add(complete_model) <- "\n"
    if (is.character(save_model)) {
      cat(complete_model, file = save_model)
    }
  }
  class(complete_model) <- c("character", "brmsmodel")
  complete_model
}
