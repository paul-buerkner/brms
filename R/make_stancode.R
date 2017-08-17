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
                          prior = NULL, autocor = NULL, nonlinear = NULL,
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
  formula <- amend_formula(
    formula, data = data, family = family, 
    autocor = autocor, threshold = threshold, 
    nonlinear = nonlinear
  )
  bterms <- parse_bf(formula)
  family <- bterms$family
  autocor <- bterms$autocor
  sample_prior <- check_sample_prior(sample_prior)
  prior <- check_prior(
    prior, formula, data = data, sample_prior = sample_prior,
    warn = !isTRUE(dots$brm_call)
  )
  data <- update_data(data, bterms = bterms)
  
  # flags to indicate the family type
  is_categorical <- is_categorical(family)
  is_mv <- is_linear(family) && length(bterms$response) > 1L
  is_forked <- is_forked(family)
  bounds <- get_bounds(bterms$adforms$trunc, data = data)
  
  ranef <- tidy_ranef(bterms, data = data)
  if (length(bterms$response) > 1L) {
    text_effects_mv <- stan_effects_mv(
      bterms, data = data, ranef = ranef, 
      prior = prior, sparse = sparse
    )
    bterms$dpars[["mu"]] <- NULL
  } else {
    text_effects_mv <- list()
  }
  text_effects <- stan_effects(
    bterms, data = data, ranef = ranef, 
    prior = prior, sparse = sparse
  )
  text_effects <- collapse_lists(text_effects, text_effects_mv)
  # because of the ID syntax, group-level effects are evaluated separately
  text_ranef <- collapse_lists(ls = lapply(
    X = unique(ranef$id), FUN = stan_re,
    ranef = ranef, prior = prior, cov_ranef = cov_ranef
  ))
  
  # generate Stan code of the likelihood
  text_llh <- stan_llh(bterms, data = data)
  # generate Stan code specific to certain models
  text_autocor <- stan_autocor(bterms, prior = prior)
  text_mv <- stan_mv(bterms, prior = prior)
  text_ordinal <- stan_ordinal(bterms, prior = prior)
  text_families <- stan_families(bterms)
  text_mixture <- stan_mixture(bterms, prior = prior)
  text_se <- stan_se(bterms)
  text_cens <- stan_cens(bterms, data)
  text_disp <- stan_disp(bterms)
  kron <- stan_needs_kronecker(ranef, names_cov_ranef = names(cov_ranef))
  text_misc_funs <- stan_misc_functions(bterms, prior, kronecker = kron)
  text_pred_funs <- stan_pred_functions(text_effects)
    
  # get priors for all parameters in the model
  text_prior <- paste0(
    text_effects$prior,
    text_ranef$prior,
    text_ordinal$prior,
    text_autocor$prior,
    text_mv$prior,
    text_mixture$prior,
    stan_prior(class = "", prior = prior)
  )
  
  # generate functions block
  text_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions { \n",
      text_misc_funs,
      text_pred_funs,
      text_autocor$fun,
      text_ordinal$fun,
      text_families$fun,
      stan_funs,
    "} \n"
  )
  
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
    text_effects$data,
    text_ranef$data,
    text_ordinal$data,
    text_families$data,
    text_mixture$data,
    text_autocor$data,
    text_cens$data,
    text_disp$data,
    text_se$data,
    if (has_trials(family))
      "  int trials[N];  // number of trials \n",
    if (is.formula(bterms$adforms$weights))
      "  vector<lower=0>[N] weights;  // model weights \n",
    if (is.formula(bterms$adforms$dec))
      "  int<lower=0,upper=1> dec[N];  // decisions \n",
    if (any(bounds$lb > -Inf))
      paste0("  ", rtype, " lb[N];  // lower truncation bounds; \n"),
    if (any(bounds$ub < Inf))
      paste0("  ", rtype, " ub[N];  // upper truncation bounds \n"),
    "  int prior_only;  // should the likelihood be ignored? \n",
    "} \n"
  )
  
  # generate transformed parameters block
  text_transformed_data <- paste0(
    "transformed data { \n",
       text_families$tdataD,
       text_effects$tdataD,
       text_se$tdataD,
       text_autocor$tdataD,
       text_families$tdataC,
       text_effects$tdataC,
       text_se$tdataC,
       text_autocor$tdataC,
    "} \n"
  )
  
  # generate parameters block
  text_parameters <- paste0(
    text_effects$par,
    text_ranef$par,
    text_ordinal$par,
    text_autocor$par,
    text_mv$par,
    text_mixture$par
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
      text_mixture$tparD,
      text_ranef$tparD,
      text_autocor$tparD, 
      text_ordinal$tparD,
      text_mv$tparD,
      text_effects$tparC1,
      text_mixture$tparC1,
      text_ranef$tparC1,
      text_autocor$tparC1, 
      text_ordinal$tparC1, 
      text_mv$tparC1,
    "} \n"
  )
  
  # generate model block
  # list dpars before pred as part of fixing issue #124
  text_model_loop <- paste0(
    text_effects$modelC2, 
    text_autocor$modelC2,
    text_effects$modelC3,
    text_mixture$modelC3,
    text_effects$modelC4
  )
  if (isTRUE(nzchar(text_model_loop))) {
    text_model_loop <- paste0(
      "  for (n in 1:N) { \n", text_model_loop, "  } \n"
    )
  }
  text_lp_pre <- list()
  if (is.formula(bterms$adforms$weights) && 
      !is.formula(bterms$adforms$cens)) {
    text_lp_pre <- list(
      modelD = "  vector[N] lp_pre; \n",
      modelC = "    target += dot_product(weights, lp_pre); \n"
    )
  }
  text_model <- paste0(
    "model { \n",
      text_effects$modelD,
      text_mixture$modelD,
      text_disp$modelD,
      text_autocor$modelD,
      text_families$modelD,
      text_lp_pre$modelD,
      text_effects$modelC1,
      text_effects$modelCgp1,
      text_mixture$modelC1,
      text_autocor$modelC1, 
      text_disp$modelC1,
      text_model_loop,
      text_families$modelC,
      "  // priors including all constants \n", 
      text_prior, 
      "  // likelihood including all constants \n",
      "  if (!prior_only) { \n  ",
      text_llh, 
      text_lp_pre$modelC,
      "  } \n", 
      text_rngprior$model,
    "} \n"
  )
  
  # generate generated quantities block
  text_generated_quantities <- paste0(
    "generated quantities { \n",
      text_effects$genD,
      text_ranef$genD,
      text_mv$genD, 
      text_rngprior$genD,
      text_effects$genC,
      text_ranef$genC,
      text_mv$genC, 
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
