
pull <- function(x, names) x[names]


# build_new_model
#
build_new_model <- function(call){

  # assign('dbg_call', call, .GlobalEnv)
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
  # assign('dbg_compile_args', compile_args, .GlobalEnv)
  model <- do_call(compile_model, compile_args)

  return(nlist(x, sdata, backend = call$backend, model, exclude, x_from_file, needs_refit))

}


# build_new_model
#
build_new_model_backup<- function(brm_call_list){

  # Extract parameters
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
  init       <- brm_call_list$init
  stan_model_args <-  brm_call_list$stan_model_args
  empty      <-  brm_call_list$empty
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
  # file_compress <- brm_call_list$file_compress
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
