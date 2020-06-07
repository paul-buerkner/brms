#' Stan Code for \pkg{brms} Models
#' 
#' Generate Stan code for \pkg{brms} models
#' 
#' @inheritParams brm
#' @param ... Other arguments for internal usage only.
#' 
#' @return A character string containing the fully commented \pkg{Stan} code 
#'   to fit a \pkg{brms} model.
#'  
#' @examples 
#' make_stancode(rating ~ treat + period + carry + (1|subject), 
#'               data = inhaler, family = "cumulative")
#' 
#' make_stancode(count ~ zAge + zBase * Trt + (1|patient),
#'               data = epilepsy, family = "poisson")
#'
#' @export
make_stancode <- function(formula, data, family = gaussian(), 
                          prior = NULL, autocor = NULL,
                          cov_ranef = NULL, sparse = NULL, 
                          sample_prior = "no", 
                          stanvars = NULL, stan_funs = NULL, 
                          knots = NULL, save_model = NULL, 
                          ...) {
  
  if (is.brmsfit(formula)) {
    stop2("Use 'stancode' to extract Stan code from 'brmsfit' objects.")
  }
  formula <- validate_formula(
    formula, data = data, family = family, 
    autocor = autocor, sparse = sparse,
    cov_ranef = cov_ranef
  )
  bterms <- brmsterms(formula)
  data <- validate_data(data, bterms = bterms, knots = knots)
  prior <- validate_prior(
    prior, bterms = bterms, data = data,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars, stan_funs = stan_funs)
  
 .make_stancode(
   bterms, data = data, prior = prior, 
   stanvars = stanvars, save_model = save_model,
   ...
 ) 
}

# internal work function of 'make_stancode'
# @param parse parse the Stan model for automatic syntax checking
# @param backend name of the backend used for parsing
# @param silent silence parsing messages
.make_stancode <- function(bterms, data, prior, stanvars,
                           parse = getOption("parse_stancode", FALSE), 
                           backend = getOption("stan_backend", "rstan"),
                           silent = TRUE, save_model = NULL, ...) {
 
  parse <- as_one_logical(parse)
  backend <- match.arg(backend, backend_choices())
  silent <- as_one_logical(silent)
  ranef <- tidy_ranef(bterms, data = data)
  meef <- tidy_meef(bterms, data = data)
  scode_predictor <- stan_predictor(
    bterms, data = data, prior = prior, 
    ranef = ranef, meef = meef,
    stanvars = stanvars
  )
  scode_ranef <- stan_re(ranef, prior = prior)
  scode_llh <- stan_llh(bterms, data = data)
  scode_global_defs <- stan_global_defs(bterms, prior = prior, ranef = ranef)
  scode_Xme <- stan_Xme(meef, prior = prior)
    
  # get priors for all parameters in the model
  scode_prior <- paste0(
    scode_predictor$prior,
    scode_ranef$prior,
    scode_Xme$prior,
    stan_unchecked_prior(prior)
  )
  # generate functions block
  scode_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions {\n",
      scode_global_defs$fun,
      collapse_stanvars(stanvars, "functions"),
    "}\n"
  )
  # generate data block
  scode_data <- paste0(
    "data {\n",
    scode_predictor$data,
    scode_ranef$data,
    scode_Xme$data,
    "  int prior_only;  // should the likelihood be ignored?\n",
    collapse_stanvars(stanvars, "data"),
    "}\n"
  )
  # generate transformed parameters block
  scode_transformed_data <- paste0(
    "transformed data {\n",
       scode_global_defs$tdata_def,
       scode_predictor$tdata_def,
       collapse_stanvars(stanvars, "tdata", "start"),
       scode_predictor$tdata_comp,
       collapse_stanvars(stanvars, "tdata", "end"),
    "}\n"
  )
  # generate parameters block
  scode_parameters <- paste0(
    scode_predictor$par,
    scode_ranef$par,
    scode_Xme$par
  )
  scode_rngprior <- stan_rngprior(
    prior = scode_prior,
    par_declars = scode_parameters,
    gen_quantities = scode_predictor$gen_def,
    prior_special = attr(prior, "special"),
    sample_prior = get_sample_prior(prior)
  )
  scode_parameters <- paste0(
    "parameters {\n",
      scode_parameters,
      scode_rngprior$par,
      collapse_stanvars(stanvars, "parameters"),
    "}\n"
  )
  # generate transformed parameters block
  scode_transformed_parameters <- paste0(
    "transformed parameters {\n",
      scode_predictor$tpar_def,
      scode_ranef$tpar_def,
      scode_Xme$tpar_def,
      collapse_stanvars(stanvars, "tparameters", "start"),
      scode_predictor$tpar_prior,
      scode_ranef$tpar_prior,
      scode_Xme$tpar_prior,
      scode_predictor$tpar_comp,
      scode_ranef$tpar_comp,
      scode_Xme$tpar_comp,
      collapse_stanvars(stanvars, "tparameters", "end"),
    "}\n"
  )
  # generate model block
  scode_model <- paste0(
    "model {\n",
      scode_predictor$model_def,
      collapse_stanvars(stanvars, "model", "start"),
      scode_predictor$model_comp_basic,
      scode_predictor$model_comp_eta_loop,
      scode_predictor$model_comp_dpar_link,
      scode_predictor$model_comp_mu_link,
      scode_predictor$model_comp_dpar_trans,
      scode_predictor$model_comp_mix,
      scode_predictor$model_comp_arma,
      scode_predictor$model_comp_catjoin,
      scode_predictor$model_comp_mvjoin,
      "  // priors including all constants\n", 
      scode_prior, 
      "  // likelihood including all constants\n",
      "  if (!prior_only) {\n",
      scode_llh, 
      "  }\n", 
      scode_rngprior$model,
      collapse_stanvars(stanvars, "model", "end"),
    "}\n"
  )
  # generate generated quantities block
  scode_generated_quantities <- paste0(
    "generated quantities {\n",
      scode_predictor$gen_def,
      scode_ranef$gen_def,
      scode_Xme$gen_def,
      scode_rngprior$gen_def,
      collapse_stanvars(stanvars, "genquant", "start"),
      scode_predictor$gen_comp,
      scode_ranef$gen_comp,
      scode_rngprior$gen_comp,
      scode_Xme$gen_comp,
      collapse_stanvars(stanvars, "genquant", "end"),
    "}\n"
  )
  # combine all elements into a complete Stan model
  scode <- paste0(
    scode_functions,
    scode_data, 
    scode_transformed_data, 
    scode_parameters,
    scode_transformed_parameters,
    scode_model,
    scode_generated_quantities
  )
  
  scode <- expand_include_statements(scode)
  if (parse) {
    scode <- parse_model(scode, backend, silent = silent)
  }
  if (is.character(save_model)) {
    cat(scode, file = save_model)
  }
  class(scode) <- c("character", "brmsmodel")
  scode
}

#' @export
print.brmsmodel <- function(x, ...) {
  cat(x)
  invisible(x) 
}

#' Extract Stan model code
#' 
#' Extract Stan code that was used to specify the model.
#' 
#' @aliases stancode.brmsfit
#' 
#' @param object An object of class \code{brmsfit}.
#' @param version Logical; indicates if the first line containing
#'   the \pkg{brms} version number should be included.
#'   Defaults to \code{TRUE}.
#' @param ... Currently ignored.
#' 
#' @return Stan model code for further processing.
#' 
#' @export
stancode.brmsfit <- function(object, version = TRUE, ...) {
  out <- object$model
  if (!version) {
    out <- sub("^[^\n]+[[:digit:]]\\.[^\n]+\n", "", out) 
  }
  out
}

#' @rdname stancode.brmsfit
#' @export
stancode <- function(object, ...) {
  UseMethod("stancode")
}

# expand '#include' statements
# This could also be done automatically by Stan at compilation time
# but would result in Stan code that is not self-contained until compilation
# @param model Stan code potentially including '#include' statements
# @return Stan code with '#include' statements expanded
expand_include_statements <- function(model) {
  path <- system.file("chunks", package = "brms")
  includes <- get_matches("#include '[^']+'", model)
  # removal of duplicates could make code generation easier in the future
  includes <- unique(includes)
  files <- gsub("(#include )|(')", "", includes)
  for (i in seq_along(includes)) {
    code <- readLines(paste0(path, "/", files[i]))
    code <- paste0(code, collapse = "\n")
    pattern <- paste0(" *", escape_all(includes[i]))
    model <- sub(pattern, code, model)
  }
  model
}
