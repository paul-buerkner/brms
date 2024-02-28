#' @title Stan Code for a Bayesian multilevel model
#'
#' @description \code{make_stancode} is a generic function that can be used to
#'   generate the Stan code to fit Bayesian multilevel models from various packages
#'   with Stan. The function invokes particular methods which depend on the
#'   class of the formula object.
#'
#'   You can view the available methods by typing:
#'
#'   \code{methods(make_stancode)}
#'
#'   See \code{\link[brms:make_stancode.default]{make_stancode}} for the default
#'   method applied for \pkg{brms}
#'
#' @param formula A formula object whose class will determine which method will
#'   be used. A symbolic description of the model to be fitted.
#' @param data An object of class data.frame, or one that can be coerced to that
#'   class) containing data of all variables used in the model.
#' @param ... Further arguments passed to the specific method
#'
#' @return A named list of objects containing the required data to fit a
#'   \pkg{brms} model with \pkg{Stan}.
#'
#' @examples
#' make_stancode(rating ~ treat + period + carry + (1|subject),
#'               data = inhaler, family = "cumulative")
#'
#' ## for more examples, see ?make_stancode.default for \pkg{brms} and for the other
#' ## methods by first calling:
#' methods(make_stancode)
#'
#' ## and then ?make_stancode.* where * is the method name
#'
#' @seealso \code{\link{make_stancode.default}}
#' @export
make_stancode <- function(formula, data, ...) {
  UseMethod("make_stancode")
}



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
make_stancode.default <- function(formula, data, family = gaussian(),
                                  prior = NULL, autocor = NULL, data2 = NULL,
                                  cov_ranef = NULL, sparse = NULL,
                                  sample_prior = "no", stanvars = NULL,
                                  stan_funs = NULL, knots = NULL,
                                  drop_unused_levels = TRUE,
                                  threads = getOption("brms.threads", NULL),
                                  normalize = getOption("brms.normalize", TRUE),
                                  save_model = NULL, ...) {

  if (is.brmsfit(formula)) {
    stop2("Use 'stancode' to extract Stan code from 'brmsfit' objects.")
  }
  formula <- validate_formula(
    formula, data = data, family = family,
    autocor = autocor, sparse = sparse,
    cov_ranef = cov_ranef
  )
  bterms <- brmsterms(formula)
  data2 <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor(formula),
    get_data2_cov_ranef(formula)
  )
  data <- validate_data(
    data, bterms = bterms,
    data2 = data2, knots = knots,
    drop_unused_levels = drop_unused_levels
  )
  prior <- .validate_prior(
    prior, bterms = bterms, data = data,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars, stan_funs = stan_funs)
  threads <- validate_threads(threads)

 .make_stancode(
   bterms, data = data, prior = prior,
   stanvars = stanvars, threads = threads,
   normalize = normalize, save_model = save_model,
   ...
 )
}

# internal work function of 'make_stancode'
# @param parse parse the Stan model for automatic syntax checking
# @param backend name of the backend used for parsing
# @param silent silence parsing messages
.make_stancode <- function(bterms, data, prior, stanvars,
                           threads = threading(),
                           normalize = getOption("brms.normalize", TRUE),
                           parse = getOption("brms.parse_stancode", FALSE),
                           backend = getOption("brms.backend", "rstan"),
                           silent = TRUE, save_model = NULL, ...) {
  normalize <- as_one_logical(normalize)
  parse <- as_one_logical(parse)
  backend <- match.arg(backend, backend_choices())
  silent <- as_one_logical(silent)
  ranef <- tidy_ranef(bterms, data = data)
  meef <- tidy_meef(bterms, data = data)
  scode_predictor <- stan_predictor(
    bterms, data = data, prior = prior,
    normalize = normalize, ranef = ranef, meef = meef,
    stanvars = stanvars, threads = threads
  )
  scode_ranef <- stan_re(
    ranef, prior = prior, threads = threads, normalize = normalize
  )
  scode_Xme <- stan_Xme(
    meef, prior = prior, threads = threads, normalize = normalize
  )
  scode_global_defs <- stan_global_defs(
    bterms, prior = prior, ranef = ranef, threads = threads
  )

  # extend Stan's likelihood part
  if (use_threading(threads)) {
    # threading is activated
    for (i in seq_along(scode_predictor)) {
      resp <- usc(names(scode_predictor)[i])
      pll_args <- stan_clean_pll_args(
        scode_predictor[[i]][["pll_args"]],
        scode_ranef[["pll_args"]],
        scode_Xme[["pll_args"]],
        collapse_stanvars_pll_args(stanvars)
      )
      partial_log_lik <- paste0(
        scode_predictor[[i]][["pll_def"]],
        scode_predictor[[i]][["model_def"]],
        collapse_stanvars(stanvars, "likelihood", "start"),
        scode_predictor[[i]][["model_comp_basic"]],
        scode_predictor[[i]][["model_comp_eta"]],
        scode_predictor[[i]][["model_comp_eta_loop"]],
        scode_predictor[[i]][["model_comp_dpar_link"]],
        scode_predictor[[i]][["model_comp_dpar_trans"]],
        scode_predictor[[i]][["model_comp_mix"]],
        scode_predictor[[i]][["model_comp_arma"]],
        scode_predictor[[i]][["model_comp_catjoin"]],
        scode_predictor[[i]][["model_comp_mvjoin"]],
        scode_predictor[[i]][["model_log_lik"]],
        collapse_stanvars(stanvars, "likelihood", "end")
      )
      partial_log_lik <- gsub(" target \\+=", " ptarget +=", partial_log_lik)
      partial_log_lik <- paste0(
        "// compute partial sums of the log-likelihood\n",
        "real partial_log_lik", resp, "_lpmf(array[] int seq", resp,
        ", int start, int end", pll_args$typed, ") {\n",
        "  real ptarget = 0;\n",
        "  int N = end - start + 1;\n",
        partial_log_lik,
        "  return ptarget;\n",
        "}\n"
      )
      partial_log_lik <- wsp_per_line(partial_log_lik, 2)
      scode_predictor[[i]][["partial_log_lik"]] <- partial_log_lik
      static <- str_if(threads$static, "_static")
      scode_predictor[[i]][["model_lik"]] <- paste0(
        "  target += reduce_sum", static, "(partial_log_lik", resp, "_lpmf",
        ", seq", resp, ", grainsize", pll_args$plain, ");\n"
      )
      str_add(scode_predictor[[i]][["tdata_def"]]) <- glue(
        "  array[N{resp}] int seq{resp} = sequence(1, N{resp});\n"
      )
    }
    scode_predictor <- collapse_lists(ls = scode_predictor)
    scode_predictor[["model_lik"]] <- paste0(
      scode_predictor[["model_no_pll_def"]],
      scode_predictor[["model_no_pll_comp_basic"]],
      scode_predictor[["model_no_pll_comp_mvjoin"]],
      scode_predictor[["model_lik"]]
    )
    str_add(scode_predictor[["data"]]) <-
      "  int grainsize;  // grainsize for threading\n"
  } else {
    # threading is not activated
    scode_predictor <- collapse_lists(ls = scode_predictor)
    scode_predictor[["model_lik"]] <- paste0(
      scode_predictor[["model_no_pll_def"]],
      scode_predictor[["model_def"]],
      collapse_stanvars(stanvars, "likelihood", "start"),
      scode_predictor[["model_no_pll_comp_basic"]],
      scode_predictor[["model_comp_basic"]],
      scode_predictor[["model_comp_eta"]],
      scode_predictor[["model_comp_eta_loop"]],
      scode_predictor[["model_comp_dpar_link"]],
      scode_predictor[["model_comp_dpar_trans"]],
      scode_predictor[["model_comp_mix"]],
      scode_predictor[["model_comp_arma"]],
      scode_predictor[["model_comp_catjoin"]],
      scode_predictor[["model_no_pll_comp_mvjoin"]],
      scode_predictor[["model_comp_mvjoin"]],
      scode_predictor[["model_log_lik"]],
      collapse_stanvars(stanvars, "likelihood", "end")
    )
  }
  scode_predictor[["model_lik"]] <-
    wsp_per_line(scode_predictor[["model_lik"]], 2)

  # get all priors added to 'lprior'
  scode_tpar_prior <- paste0(
    scode_predictor[["tpar_prior"]],
    scode_ranef[["tpar_prior"]],
    scode_Xme[["tpar_prior"]]
  )

  # generate functions block
  scode_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions {\n",
      scode_global_defs[["fun"]],
      collapse_stanvars(stanvars, "functions"),
      scode_predictor[["partial_log_lik"]],
    "}\n"
  )

  # generate data block
  scode_data <- paste0(
    "data {\n",
    "  int<lower=1> N;  // total number of observations\n",
    scode_predictor[["data"]],
    scode_ranef[["data"]],
    scode_Xme[["data"]],
    "  int prior_only;  // should the likelihood be ignored?\n",
    collapse_stanvars(stanvars, "data"),
    "}\n"
  )

  # generate transformed parameters block
  scode_transformed_data <- paste0(
    "transformed data {\n",
       scode_global_defs[["tdata_def"]],
       scode_predictor[["tdata_def"]],
       collapse_stanvars(stanvars, "tdata", "start"),
       scode_predictor[["tdata_comp"]],
       collapse_stanvars(stanvars, "tdata", "end"),
    "}\n"
  )

  # generate parameters block
  scode_parameters <- paste0(
    scode_predictor[["par"]],
    scode_ranef[["par"]],
    scode_Xme[["par"]]
  )
  # prepare additional sampling from priors
  scode_rngprior <- stan_rngprior(
    tpar_prior = scode_tpar_prior,
    par_declars = scode_parameters,
    gen_quantities = scode_predictor[["gen_def"]],
    special_prior = attr(prior, "special"),
    sample_prior = get_sample_prior(prior)
  )
  scode_parameters <- paste0(
    "parameters {\n",
      scode_parameters,
      scode_rngprior[["par"]],
      collapse_stanvars(stanvars, "parameters"),
    "}\n"
  )

  # generate transformed parameters block
  scode_lprior_def <- "  real lprior = 0;  // prior contributions to the log posterior\n"
  scode_transformed_parameters <- paste0(
    "transformed parameters {\n",
      scode_predictor[["tpar_def"]],
      scode_ranef[["tpar_def"]],
      scode_Xme[["tpar_def"]],
      str_if(normalize, scode_lprior_def),
      collapse_stanvars(stanvars, "tparameters", "start"),
      scode_predictor[["tpar_prior_const"]],
      scode_ranef[["tpar_prior_const"]],
      scode_Xme[["tpar_prior_const"]],
      scode_predictor[["tpar_comp"]],
      scode_predictor[["tpar_special_prior"]],
      scode_ranef[["tpar_comp"]],
      scode_Xme[["tpar_comp"]],
      # lprior cannot contain _lupdf functions in transformed parameters
      # as discussed on github.com/stan-dev/stan/issues/3094
      str_if(normalize, scode_tpar_prior),
      collapse_stanvars(stanvars, "tparameters", "end"),
    "}\n"
  )

  # combine likelihood with prior part
  not_const <- str_if(!normalize, " not")
  scode_model <- paste0(
    "model {\n",
      str_if(!normalize, scode_lprior_def),
      collapse_stanvars(stanvars, "model", "start"),
      "  // likelihood", not_const, " including constants\n",
      "  if (!prior_only) {\n",
      scode_predictor[["model_lik"]],
      "  }\n",
      "  // priors", not_const, " including constants\n",
      str_if(!normalize, scode_tpar_prior),
      "  target += lprior;\n",
      scode_predictor[["model_prior"]],
      scode_ranef[["model_prior"]],
      scode_Xme[["model_prior"]],
      stan_unchecked_prior(prior),
      collapse_stanvars(stanvars, "model", "end"),
    "}\n"
  )
  # generate generated quantities block
  scode_generated_quantities <- paste0(
    "generated quantities {\n",
      scode_predictor[["gen_def"]],
      scode_ranef[["gen_def"]],
      scode_Xme[["gen_def"]],
      scode_rngprior[["gen_def"]],
      collapse_stanvars(stanvars, "genquant", "start"),
      scode_predictor[["gen_comp"]],
      scode_ranef[["gen_comp"]],
      scode_rngprior[["gen_comp"]],
      scode_Xme[["gen_comp"]],
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
  # if (backend == "cmdstanr") {
  #   if (requireNamespace("cmdstanr", quietly = TRUE) &&
  #       cmdstanr::cmdstan_version() >= "2.29.0") {
  #     tmp_file <- cmdstanr::write_stan_file(scode)
  #     scode <- .canonicalize_stan_model(tmp_file, overwrite_file = FALSE)
  #   }
  # }
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
#' @param version Logical; indicates if the first line containing the \pkg{brms}
#'   version number should be included. Defaults to \code{TRUE}.
#' @param regenerate Logical; indicates if the Stan code should be regenerated
#'   with the current \pkg{brms} version. By default, \code{regenerate} will be
#'   \code{FALSE} unless required to be \code{TRUE} by other arguments.
#' @param threads Controls whether the Stan code should be threaded. See
#'   \code{\link{threading}} for details.
#' @param backend Controls the Stan backend. See \code{\link{brm}} for details.
#' @param ... Further arguments passed to
#'   \code{\link[brms:make_stancode.default]{make_stancode}} if the Stan code is
#'   regenerated.
#'
#' @return Stan model code for further processing.
#'
#' @export
stancode.brmsfit <- function(object, version = TRUE, regenerate = NULL,
                             threads = NULL, backend = NULL, ...) {
  if (is.null(regenerate)) {
    # determine whether regenerating the Stan code is required
    regenerate <- FALSE
    cl <- match.call()
    if ("threads" %in% names(cl)) {
      threads <- validate_threads(threads)
      if (use_threading(threads) && !use_threading(object$threads) ||
          !use_threading(threads) && use_threading(object$threads)) {
        # threading changed; regenerated Stan code
        regenerate <- TRUE
      }
      object$threads <- threads
    }
    if ("backend" %in% names(cl)) {
      backend <- match.arg(backend, backend_choices())
      # older Stan versions do not support array syntax
      if (require_old_stan_syntax(object, backend, "2.29.0")) {
        regenerate <- TRUE
      }
      object$backend <- backend
    }
  }
  regenerate <- as_one_logical(regenerate)
  if (regenerate) {
    object <- restructure(object)
    out <- make_stancode(
      formula = object$formula,
      data = object$data,
      prior = object$prior,
      data2 = object$data2,
      stanvars = object$stanvars,
      sample_prior = get_sample_prior(object$prior),
      threads = object$threads,
      backend = object$backend,
      ...
    )
  } else {
    # extract Stan code unaltered
    out <- object$model
  }
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

# check if Stan code includes normalization constants
is_normalized <- function(stancode) {
  !grepl("_lup(d|m)f\\(", stancode)
}

# Normalizes Stan code to avoid triggering refit after whitespace and
# comment changes in the generated code.
# In some distant future, StanC3 may provide its own normalizing functions,
# until then this is a set of regex hacks.
# @param x a string containing the Stan code
normalize_stancode <- function(x) {
  x <- as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+"," ", x)
  trimws(x)
}

# check if the currently installed Stan version requires older syntax
# than the Stan version with which the model was initially fitted
require_old_stan_syntax <- function(object, backend, version) {
  stopifnot(is.brmsfit(object))
  isTRUE(
    (object$backend == "rstan" && object$version$rstan >= version ||
       object$backend == "cmdstanr" && object$version$cmdstan >= version) &&
      (backend == "rstan" && utils::packageVersion("rstan") < version ||
         backend == "cmdstanr" && cmdstanr::cmdstan_version() < version)
  )
}
