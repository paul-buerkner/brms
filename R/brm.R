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

#' Validate a `brm_call` object
#'
#' Internal helper – stops with an informative error if core fields
#' are missing or malformed.  Returns the input invisibly on success.
#'
#' Uses **rlang** predicates for concise, CRAN-friendly checks.
#'
#' @param brm_call A list with class <brm_call>.
#' @return `brm_call` (invisibly) or an error via `rlang::abort()`.
#' @noRd
brm_call_type_check <- function(brm_call) {

  # --------------------------------------------------------------------- #
  # 1.  Top-level class ---------------------------------------------------
  if (!inherits(brm_call, "brm_call")) {
    rlang::abort("`brm_call` must inherit from class <brm_call>.",
                 arg = "brm_call")
  }

  # --------------------------------------------------------------------- #
  # 2.  Required fields  --------------------------------------------------
  req <- list(
    formula = function(x) is.formula(x) ||  is.brmsformula(x) || is.list(x),

    # TODO Currently we let data to be null.
    # We may modify this part later. Some of the parameters not checked yet.
    # data    = is.data.frame,
    # family  = function(x) is.family(x) || rlang::is_string(x),
    backend = rlang::is_string,
    iter    = rlang::is_scalar_integerish,
    chains  = rlang::is_scalar_integerish,
    call_only = rlang::is_scalar_logical
  )

  for (nm in names(req)) {
    if (!nm %in% names(brm_call)) {
      rlang::abort(
        glue::glue("Field `{nm}` is missing from `brm_call`."),
        arg = nm
      )
    }
    if (!isTRUE(req[[nm]](brm_call[[nm]]))) {
      rlang::abort(
        glue::glue("Field `{nm}` has the wrong type or length."),
        arg = nm
      )
    }
  }

  # --------------------------------------------------------------------- #
  # 3.  Value constraints  -----------------------------------------------
  if (brm_call$iter <= 0) {
    rlang::abort("`iter` must be a positive integer.", arg = "iter")
  }
  if (brm_call$chains <= 0) {
    rlang::abort("`chains` must be a positive integer.", arg = "chains")
  }
  if (rlang::is_true(brm_call$call_only) && !is.null(brm_call$fit)) {
    rlang::abort("`brm_call` with `call_only = TRUE` must not contain a `fit`.",
                 arg = "fit")
  }

  invisible(brm_call)
}


#' Internal engine to evaluate and fit a brms model
#' @noRd
.brm_internal <- function(brm_call) {

  # Check if fit object can be reused from file
  result <- .brm_check(brm_call)
  if (!is.null(result)) {
    return(result)
  }

  # initialize brmsfit object
  if (is.brmsfit(brm_call$fit)) {
    # re-use existing model
    .list    <- re_use_existing_model(brm_call)
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
    .list   <-  build_new_model(brm_call)
    backend <- .list$backend
    model   <- .list$model
    exclude <- .list$exclude
    x       <- .list$x
    sdata   <- .list$sdata
  }

  if(brm_call$empty){
    return(.list$x)
  }
  # ==================================
  # model, exclude, backend, x, sdata may be changed or created
  # ===================================
  fit_args <- c(
    nlist(
      model, sdata, # maybe modifed above
      algorithm = brm_call$algorithm,
      backend, # maybe modifed above
      iter = brm_call$iter, warmup = brm_call$warmup, thin = brm_call$thin,
      chains = brm_call$chains, cores = brm_call$cores,
      threads = brm_call$threads, opencl = brm_call$opencl,
      init = brm_call$init,
      exclude, # maybe modifed above
      control = brm_call$control, future = brm_call$future,
      seed = brm_call$seed, silent = brm_call$silent
    ),
    brm_call$dot_args
  )

  x$fit <- do_call(fit_model, fit_args)

  # rename parameters to have human readable names
  if (brm_call$rename) {
    x <- rename_pars(x)
  }
  if (!is.null(brm_call$file)) {
    x <- write_brmsfit(x, brm_call$file, compress = brm_call$file_compress)
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
                empty = FALSE, rename = TRUE, call_only = FALSE, ...) {

  # if called with a `brm_call` object handle it first
  if(inherits(formula , 'brm_call')){
    args <- formula
    return( .brm_internal(args))
  }

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
  call_only <- as_one_logical(call_only)
  args <- .brm_collect_args(...)
  class(args) <- c("brm_call" , "list")
  # if call_only is TRUE return a brm_call object
  if(call_only){
    return(args)
  }

  brm_call_type_check(args)
  .brm_internal(args)
}
