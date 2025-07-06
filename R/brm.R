#' Check for an existing cached `brmsfit` and return it if valid
#'
#' Early-exit helper used by `.brm_internal()`.
#' If the user supplied a `file` argument and the cached fit can be reused
#' under the chosen `file_refit` policy, this function loads the object and
#' hands it back; otherwise it returns `NULL` and the caller proceeds to
#' (re-)build the model.
#'
#' **`file_refit` rules**
#' * `"never"`   – always reuse the fit (default).
#' * `"always"`  – never reuse; force a refit.
#' * `"on_change"` – reuse only when the *hash* of the current call matches
#'   the hash stored in the cached object (attribute `"brm_call_hash"`).
#'
#' @param call A validated **`brm_call`** list.
#' @return A `brmsfit` object **or** `NULL` if no valid cache can be used.
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

#' Re-use a compiled Stan model when possible
#'
#' Internal helper called from `.brm_internal()` when the user passed a
#' **previously fitted** `brmsfit` object via the `fit` field of the
#' `brm_call`.  It decides whether that compiled model can be recycled,
#' pulls out the pieces we still need, and returns them in a tidy bundle.
#'
#' @param call A validated **`brm_call`** list whose `fit` element is a
#'   `brmsfit` object.
#'
#' @return A named list with elements
#'   \describe{
#'     \item{`backend`}{Backend string (`"rstan"`, `"cmdstanr"`, or `"mock"`).}
#'     \item{`model`}{Compiled Stan model reused from the old fit.}
#'     \item{`exclude`}{Names of parameters that should be dropped after
#'       sampling (copied from the old fit).}
#'     \item{`x`}{The old `brmsfit` object, ready for updating.}
#'     \item{`sdata`}{The standata list extracted from the fit.}
#'     \item{`needs_refit`}{Logical flag; `FALSE` if the cached fit can be
#'       returned as-is, `TRUE` if we must run the sampler again.}
#'     \item{`x_from_file`}{`NULL` – kept for API symmetry with the
#'       `build_new_model()` helper.}
#'   }
#' @noRd
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

# Internal method to create fit args
#' @noRd
.create_fit_args <- function(brm_call){
  #   model, exclude, backend, x, sdata may be changed or created in `.build_or_reuse`
  .list <- .build_or_reuse(brm_call)
  backend <- .list$backend
  model   <- .list$model
  exclude <- .list$exclude
  sdata   <- .list$sdata

  # maybe modifed by .build_or_reuse
  fit_args <- c(
    nlist(
      model, sdata, # maybe modifed by .build_or_reuse
      algorithm = brm_call$algorithm,
      backend, # maybe modifed
      iter = brm_call$iter, warmup = brm_call$warmup, thin = brm_call$thin,
      chains = brm_call$chains, cores = brm_call$cores,
      threads = brm_call$threads, opencl = brm_call$opencl,
      init = brm_call$init,
      exclude, # maybe modifed
      control = brm_call$control, future = brm_call$future,
      seed = brm_call$seed, silent = brm_call$silent
    ),
    brm_call$dot_args
  )
  fit_args
}

# Internal function that will create or reuse existing brmsfit
#
.build_or_reuse <- function(brm_call){
  # ====================================================================
  #   model, exclude, backend, x, sdata may be changed or created here
  # ====================================================================
  # initialize brmsfit object
  if (is.brmsfit(brm_call$fit)) {
    # re-use existing model
    .list    <- re_use_existing_model(brm_call)
  } else {
    # build new model
    .list   <-  build_new_model(brm_call)
  }
  .list
}

#' Internal engine to evaluate and fit a *brms* model
#' @noRd
.brm_internal <- function(brm_call) {

  # Check if fit object can be reused from file
  result <- .brm_check(brm_call)
  if (!is.null(result)) {
    return(result)
  }

  # build new or reuse existing fit
  .list <- .build_or_reuse(brm_call)
  if(!.list$needs_refit){
     return(.list$x_from_file) # return x from file
  }
  # empty model returns
  if(brm_call$empty){
    return(.list$x)
  }

  # brmsfit object `x`
  x <- .list$x
  # fit happens
  fit_args <- .create_fit_args(brm_call)
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

#' Collect `brm()` arguments into a tidy **brm_call** object
#'
#' Internal helper used at the very top of `brm()`.
#' It separates *formal* `brm()` arguments from any extra
#' Stan-backend tuning options that a user might pass through `...`,
#' then stores everything in a lightweight list with class
#' **`brm_call`**.  The resulting object goes straight to
#' `brm_call_type_check()` and then into `.brm_internal()`.
#'
#' *Implementation notes*
#' * We grab the names of the **current** `brm()` formals at run-time
#'   (`names(formals(brms::brm))`) so the helper automatically stays in
#'   sync with upstream changes in **brms**.
#' * Any argument not in that set is treated as an
#'   *extra* Stan argument and saved under `dot_args`.
#'
#' @param ... Arguments passed from the public `brm()` wrapper.
#' @return A list of class `c("brm_call", "list")`.
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

  call_only <- as_one_logical(call_only)
  # if called with a `brm_call` object handle it first
  if(is.brm_call(formula)) {
    args <- formula
    if(call_only) {
      # when called with brm(brm_call, call_only = TRUE)
      # also returns `brm_call`
      return(args)
    }

    return(.brm_internal(args))
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

  args <- .brm_collect_args(...)
  class(args) <- c("brm_call" , "list")
  # if call_only is TRUE return a brm_call object
  if(call_only){
    return(args)
  }

  brm_call_type_check(args)
  .brm_internal(args)
}
