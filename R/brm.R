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
  algorithm    <- brm_call_list$algorithm
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
      model, sdata, algorithm,backend = backend, iter, warmup, thin, chains, cores,
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
  # args
  .brm_internal(args)
}
