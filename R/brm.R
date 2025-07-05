#' Check for an existing cached brmsfit and return it if valid
#' @noRd
.brm_check <- function(brm_call_list) {
  file       <- brm_call_list$file
  file_refit <- match.arg(brm_call_list$file_refit, file_refit_options())

  # Load brmsfit only if refit is explicitly set to 'never'
  if (!is.null(file) && file_refit == "never") {
    x <- read_brmsfit(file)
    if (!is.null(x)) {
      return(x)
    }
  }
  NULL
}

# re_use_existing_model
#
re_use_existing_model<- function(.brm_call_list){

  # re-use existing model
  x <- .brm_call_list$fit
  x$criteria <- list()
  sdata <- standata(x)
  x_from_file <- NULL
  if (!is.null(.brm_call_list$file) && .brm_call_list$file_refit == "on_change") {
    x_from_file <- read_brmsfit(.brm_call_list$file)
    if (!is.null(x_from_file)) {
      needs_refit <- brmsfit_needs_refit(
        x_from_file, scode = stancode(x), sdata = sdata,
        data = x$data, algorithm = .brm_call_list$algorithm, silent = .brm_call_list$silent
      )
      if (!needs_refit) {
        return(x_from_file) # TODO
      }
    }
  }
  backend <- x$backend
  model <- compiled_model(x)
  exclude <- exclude_pars(x)
  return(nlist(backend, model, exclude, x_from_file, needs_refit))
}

# build_new_model
#
build_new_model<- function(brm_call_list){

  # Extract parameters
  fit        <- brm_call_list$fit
  file       <- brm_call_list$file
  file_refit <- brm_call_list$file_refit
  algorithm  <- brm_call_list$algorithm
  backend    <- brm_call_list$backend
  silent     <- brm_call_list$silent

  formula    <- brm_call_list$formula
  data       <- brm_call_list$data
  family     <- brm_call_list$family
  autocor    <- brm_call_list$autocor
  sparse     <- brm_call_list$sparse
  cov_ranef  <- brm_call_list$cov_ranef
  data2      <- brm_call_list$data2
  knots      <- brm_call_list$knots
  drop_unused_levels <- brm_call_list$drop_unused_levels
  prior      <- brm_call_list$prior
  bframe     <- brm_call_list$bframe
  sample_prior <- brm_call_list$sample_prior
  stanvars   <- brm_call_list$stanvars
  stan_funs  <- brm_call_list$stan_funs
  save_pars  <- brm_call_list$save_pars
  save_ranef <- brm_call_list$save_ranef
  save_mevars <- brm_call_list$save_mevars
  save_all_pars <- brm_call_list$save_all_pars

  model      <- brm_call_list$model

  save_model <- brm_call_list$save_model
  backend    <- brm_call_list$backend
  threads    <- brm_call_list$threads
  opencl     <- brm_call_list$opencl

  normalize  <- brm_call_list$normalize
  control    <- brm_call_list$control



  rename     <- brm_call_list$rename
  init       <- brm_call_list$init
  stan_model_args <-  brm_call_list$stan_model_args
  empty      <-  brm_call_list$empty
  iter       <- brm_call_list$iter
  warmup     <- brm_call_list$warmup
  thin       <- brm_call_list$thin
  chains     <- brm_call_list$chains
  cores      <- brm_call_list$cores
  threads    <- brm_call_list$threads
  opencl     <- brm_call_list$opencl
  exclude    <- brm_call_list$exclude
  control    <- brm_call_list$control
  future     <- brm_call_list$future
  seed       <- brm_call_list$seed
  silent     <- brm_call_list$silent
  file_compress <- brm_call_list$file_compress

  needs_refit <- TRUE

  # --- Build a new brmsfit object from scratch ---
  # ========================================================
  # build new model
  x_from_file <- NULL
  formula <- validate_formula(
    formula, data = data, family = family,
    autocor = autocor, sparse = sparse,
    cov_ranef = cov_ranef
  )
  family <- get_element(formula, "family")
  bterms <- brmsterms(formula)
  data2  <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor(formula),
    get_data2_cov_ranef(formula)
  )
  data <- validate_data(
    data, bterms = bterms,
    data2 = data2, knots = knots,
    drop_unused_levels = drop_unused_levels,
    data_name = substitute_name(data)
  )
  bframe <- brmsframe(bterms, data)
  prior <- .validate_prior(
    prior, bframe = bframe,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars, stan_funs = stan_funs)
  save_pars <- validate_save_pars(
    save_pars, save_ranef = save_ranef,
    save_mevars = save_mevars,
    save_all_pars = save_all_pars
  )

  # generate Stan code
  model <- .stancode(
    bframe, prior = prior, stanvars = stanvars,
    save_model = save_model, backend = backend, threads = threads,
    opencl = opencl, normalize = normalize
  )
  # initialize S3 object
  stan_args <- c(
    nlist(init, silent, control, stan_model_args),
    brm_call_list$dot_args
  )
  x <- brmsfit(
    formula = formula, data = data, data2 = data2, prior = prior,
    stanvars = stanvars, model = model, algorithm = algorithm,
    backend = backend, threads = threads, opencl = opencl,
    save_pars = save_pars, ranef = bframe$frame$re, family = family,
    basis = frame_basis(bframe, data = data),
    stan_args = stan_args
  )
  exclude <- exclude_pars(x, bframe = bframe)
  # generate Stan data before compiling the model to avoid
  # unnecessary compilations in case of invalid data
  sdata <- .standata(
    bframe, data = data, prior = prior, data2 = data2,
    stanvars = stanvars, threads = threads
  )
  if (empty) {
    # return the brmsfit object with an empty 'fit' slot
    # return(x)
    model <- NULL
    x_from_file <- NULL
    return(nlist(x, sdata, backend, model, exclude, x_from_file, needs_refit))
  }

  if (!is.null(file) && file_refit == "on_change") {
    x_from_file <- read_brmsfit(file)
    if (!is.null(x_from_file)) {
      needs_refit <- brmsfit_needs_refit(
        x_from_file, scode = model, sdata = sdata, data = data,
        algorithm = algorithm, silent = silent
      )

      if (!needs_refit) {
        model <- NULL
        return(nlist(x, sdata, backend, model, exclude, x_from_file, needs_refit))
      }
    }
  }

  # compile the Stan model
  compile_args <- stan_model_args
  compile_args$model <- model
  compile_args$backend <- backend
  compile_args$threads <- threads
  compile_args$opencl <- opencl
  compile_args$silent <- silent
  model <- do_call(compile_model, compile_args)
  return(nlist(x, sdata, backend, model, exclude, x_from_file, needs_refit))

}

#' Internal engine to evaluate and fit a brms model
#' @noRd
.brm_internal <- function(brm_call_list) {

  # Extract parameters
  fit        <- brm_call_list$fit
  model      <- brm_call_list$model
  file       <- brm_call_list$file
  file_compress <- brm_call_list$file_compress
  empty        <-  brm_call_list$empty
  rename       <- brm_call_list$rename
  algorithm  <- brm_call_list$algorithm
  backend      <- brm_call_list$backend
  iter         <- brm_call_list$iter
  warmup       <- brm_call_list$warmup
  thin         <- brm_call_list$thin
  chains       <- brm_call_list$chains
  cores        <- brm_call_list$cores
  threads      <- brm_call_list$threads
  opencl       <- brm_call_list$opencl
  init         <- brm_call_list$init
  exclude      <- brm_call_list$exclude
  control      <- brm_call_list$control
  future       <- brm_call_list$future
  seed         <- brm_call_list$seed
  silent       <- brm_call_list$silent

  # Check if fit object can be reused from file
  result <- .brm_check(brm_call_list)
  if (!is.null(result)) {
    return(result)
  }

  # initialize brmsfit object
  if (is.brmsfit(fit)) {
    # re-use existing model
    .list    <- re_use_existing_model( brm_call_list )
    backend  <- .list$backend
    model    <- .list$model
    exclude  <- .list$exclude
    x        <- .list$x
    sdata    <- .list$sdata
    if(!.list$needs_refit){
      return(.list$x_from_file)
    }

  } else {
    # build new model
    .list   <-  build_new_model(brm_call_list)
    backend <- .list$backend
    model   <- .list$model
    exclude <- .list$exclude
    x       <- .list$x
    sdata   <- .list$sdata
  }

  if(empty){
    return(.list$x)
  }
  # ==================================
  # model, exclude, backend, x, sdata may be changed or created
  # ===================================
  fit_args <- c(
    nlist(
      model, sdata, algorithm, backend, iter, warmup, thin, chains, cores,
      threads, opencl, init, exclude, control, future, seed, silent
    ),
    brm_call_list$dot_args
  )

  x$fit <- do_call(fit_model, fit_args)

  # rename parameters to have human readable names
  if (rename) {
    x <- rename_pars(x)
  }
  if (!is.null(file)) {
    x <- write_brmsfit(x, file, compress = file_compress)
  }
  x
}

#' Collect brm() arguments into a clean list
#' @noRd
.brm_collect_args <- function(...) {

  call_env  <- parent.frame()
  arg_names <- names(formals(brm))

  ## 1. drop the literal "..." from arg_names
  arg_names <- arg_names[arg_names != "..."]

  ## 2. capture every formal (already evaluated inside brm())
  arg_list <- setNames(
    lapply(arg_names, function(a) get(a, envir = call_env)),
    arg_names
  )

  ## 3. stash the dot-args for later splicing
  arg_list$dot_args <- list(...)
  arg_list
}


#' @rdname brm
#' @export
brm <- function(formula, data= NULL, family = gaussian(), prior = NULL,
                autocor = NULL, data2 = NULL, cov_ranef = NULL,
                sample_prior = "no", sparse = NULL, knots = NULL,
                drop_unused_levels = TRUE, stanvars = NULL, stan_funs = NULL,
                fit = NA, save_pars = getOption("brms.save_pars", NULL),
                save_ranef = NULL, save_mevars = NULL, save_all_pars = NULL,
                init = NULL, inits = NULL, chains = 4,
                iter = getOption("brms.iter", 2000),
                warmup = floor(iter / 2), thin = 1,
                cores = getOption("mc.cores", 1),
                threads = getOption("brms.threads", NULL),
                opencl = getOption("brms.opencl", NULL),
                normalize = getOption("brms.normalize", TRUE),
                control = NULL,
                algorithm = getOption("brms.algorithm", "sampling"),
                backend = getOption("brms.backend", "rstan"),
                future = getOption("future", FALSE), silent = 1,
                seed = NA, save_model = NULL, stan_model_args = list(),
                file = NULL, file_compress = TRUE,
                file_refit = getOption("brms.file_refit", "never"),
                empty = FALSE, rename = TRUE, ...) {

  algorithm <- match.arg(algorithm, algorithm_choices())
  backend <- match.arg(backend, backend_choices())
  normalize <- as_one_logical(normalize)
  silent <- validate_silent(silent)
  iter <- as_one_numeric(iter)
  warmup <- as_one_numeric(warmup)
  thin <- as_one_numeric(thin)
  chains <- as_one_numeric(chains)
  cores <- as_one_numeric(cores)
  init <- use_alias(init, inits)
  threads <- validate_threads(threads)
  opencl <- validate_opencl(opencl)
  future <- as_one_logical(future) && chains > 0L
  seed <- as_one_numeric(seed, allow_na = TRUE)
  empty <- as_one_logical(empty)
  rename <- as_one_logical(rename)

  args <- .brm_collect_args(...)
  .brm_internal(args)
}
