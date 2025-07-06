#' Build a fresh `brmsfit` shell from a validated call
#'
#' Internal helper used by `brm()` when no existing fit can be re-used.
#' It prepares data, generates Stan code, compiles the model (unless
#' `empty = TRUE`), and returns all artefacts needed for sampling.
#'
#' @param call A list of class **`brm_call`** containing *validated* user
#'   arguments.  The object is assumed to have passed `brm_call_type_check()`.
#'
#' @return A named list with components
#'   \describe{
#'     \item{`x`}{A pre-initialised **`brmsfit`** object (may have an empty
#'       `fit` slot if `empty = TRUE`).}
#'     \item{`sdata`}{The standata list (NULL when `empty = TRUE`).}
#'     \item{`backend`}{Backend string (`"rstan"`, `"cmdstanr"`, or `"mock"`).}
#'     \item{`model`}{Compiled Stan model or `NULL` when `empty = TRUE`.}
#'     \item{`exclude`}{Character vector of parameter names to drop
#'       after sampling.}
#'     \item{`x_from_file`}{A cached fit (if `file` was supplied and no refit
#'       was needed) or `NULL`.}
#'     \item{`needs_refit`}{Logical flag used upstream to decide whether the
#'       sampler must be run.}
#'   }
#' @noRd
build_new_model <- function(call){

  needs_refit <- TRUE
  # --- Build a new brmsfit object from scratch ---
  # ========================================================
  # build new model
  x_from_file <- NULL
  formula <- validate_formula(
    call$formula, data = call$data, family = call$family,
    autocor = call$autocor, sparse = call$sparse,
    cov_ranef = call$cov_ranef
  )
  family <- get_element(formula, "family")
  bterms <- brmsterms(formula)
  data2  <- validate_data2(
    call$data2, bterms = bterms,
    get_data2_autocor(formula),
    get_data2_cov_ranef(formula)
  )

  data <- validate_data(
    call$data, bterms = bterms,
    data2 = data2, knots = call$knots,
    drop_unused_levels = call$drop_unused_levels,
    data_name = substitute_name(data)
  )

  bframe <- brmsframe(bterms, data)
  prior <- .validate_prior(
    call$prior, bframe =  bframe,
    sample_prior = call$sample_prior
  )

  stanvars <- validate_stanvars(call$stanvars, stan_funs = call$stan_funs)
  save_pars <- validate_save_pars(
    call$save_pars, save_ranef = call$save_ranef,
    save_mevars = call$save_mevars,
    save_all_pars = call$save_all_pars
  )
  # generate Stan code
  model <- .stancode(
    bframe, prior = prior, stanvars = stanvars,
    save_model = call$save_model, backend = call$backend, threads = call$threads,
    opencl = call$opencl, normalize = call$normalize
  )
  # initialize S3 object
  stan_args <- c(
    nlist(init = call$init, silent= call$silent, control = call$control, stan_model_args = call$stan_model_args),
    call$dot_args
  )
  x <- brmsfit(
    formula = formula, data = data, data2 = data2, prior = prior,
    stanvars = stanvars, model = model, algorithm = call$algorithm,
    backend = call$backend, threads = call$threads, opencl = call$opencl,
    save_pars = save_pars, ranef = bframe$frame$re, family = family,
    basis = frame_basis(bframe, data = data),
    stan_args = stan_args
  )

  x$brm_call <- call
  exclude <- exclude_pars(x, bframe = bframe)
  # generate Stan data before compiling the model to avoid
  # unnecessary compilations in case of invalid data
  sdata <- .standata(
    bframe, data = data, prior = prior, data2 = data2,
    stanvars = stanvars, threads = call$threads
  )
  if (call$empty) {
    # return the brmsfit object with an empty 'fit' slot

    model <- NULL
    x_from_file <- NULL
    return(nlist(x, sdata, call$backend, model, exclude, x_from_file, needs_refit))
  }

  if (!is.null(call$file) && call$file_refit == "on_change") {
    x_from_file <- read_brmsfit(call$file)
    if (!is.null(x_from_file)) {
      needs_refit <- brmsfit_needs_refit(
        x_from_file, scode = model, sdata = sdata, data = data,
        algorithm = call$algorithm, silent = call$silent
      )

      if (!needs_refit) {
        model <- NULL
        return(nlist(x, sdata, backend = call$backend, model, exclude, x_from_file, needs_refit))
      }
    }
  }

  # compile the Stan model
  compile_args <- c(call$stan_model_args,
                    nlist(model, backend = call$backend, threads = call$threads,
                          opencl = call$opencl, silent = call$silent))
  model <- do_call(compile_model, compile_args)

  return(nlist(x, sdata, backend = call$backend, model, exclude, x_from_file, needs_refit))
}
