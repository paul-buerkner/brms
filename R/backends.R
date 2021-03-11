# parse Stan model code
# @param model Stan model code
# @return validated Stan model code
parse_model <- function(model, backend, ...) {
  backend <- as_one_character(backend)
  .parse_model <- get(paste0(".parse_model_", backend), mode = "function")
  .parse_model(model, ...)
}

# parse Stan model code with rstan
# @param model Stan model code
# @return validated Stan model code
.parse_model_rstan <- function(model, silent = 1, ...) {
  out <- eval_silent(
    rstan::stanc(model_code = model, ...),
    type = "message", try = TRUE, silent = silent
  )
  out$model_code
}

# parse Stan model code with cmdstanr
# @param model Stan model code
# @return validated Stan model code
.parse_model_cmdstanr <- function(model, silent = 1, ...) {
  require_package("cmdstanr")
  temp_file <- cmdstanr::write_stan_file(model)
  out <- eval_silent(
    cmdstanr::cmdstan_model(temp_file, compile = FALSE, ...),
    type = "message", try = TRUE, silent = silent
  )
  out$check_syntax(quiet = TRUE)
  collapse(out$code(), "\n")
}

.parse_model_mock <- function(model, silent = TRUE, parse_error = NULL,
                              parse_check = "rstan", ...) {
  if(!is.null(parse_error)) {
    stop2(parse_error)
  } else if(parse_check == "rstan") {
    .parse_model_rstan(model, silent = silent, ...)
  } else if(parse_check == "cmdstanr") {
    .parse_model_cmdstanr(model, silent = silent, ...)
  } else if(is.null(parse_check)) {
    "mock_code"
  } else {
    stop2("Unknown parse_check value")
  }
}

# compile Stan model
# @param model Stan model code
# @return validated Stan model code
compile_model <- function(model, backend, ...) {
  backend <- as_one_character(backend)
  .compile_model <- get(paste0(".compile_model_", backend), mode = "function")
  .compile_model(model, ...)
}

# compile Stan model with rstan
# @param model Stan model code
# @return model compiled with rstan
.compile_model_rstan <- function(model, threads, silent = 1, ...) {
  args <- list(...)
  args$model_code <- model
  if (silent < 2) {
    message("Compiling Stan program...")
  }
  if (use_threading(threads)) {
    if (utils::packageVersion("rstan") >= 2.26) {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(rstan::rstan_options(threads_per_chain = threads_per_chain_def))
      rstan::rstan_options(threads_per_chain = threads$threads)
    } else {
      stop2("Threading is not supported by backend 'rstan' version ",
            utils::packageVersion("rstan"), ".")
    }
  }
  eval_silent(
    do_call(rstan::stan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

# compile Stan model with cmdstanr
# @param model Stan model code
# @return model compiled with cmdstanr
.compile_model_cmdstanr <- function(model, threads, silent = 1, ...) {
  require_package("cmdstanr")
  args <- list(...)
  args$stan_file <- cmdstanr::write_stan_file(model)
  if (use_threading(threads)) {
    args$cpp_options$stan_threads <- TRUE
  }
  eval_silent(
    do_call(cmdstanr::cmdstan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}



.compile_model_mock <- function(model, threads, compile_check = "rstan",
                                compile_error = NULL, ...) {
  if(!is.null(compile_error)) {
    stop2(compile_error)
  } else if(compile_check == "rstan") {
    .parse_model_rstan(model, silent = TRUE, ...)
  } else if(compile_check == "cmdstanr") {
    .parse_model_cmdstanr(model, silent = TRUE, ...)
  } else if(is.null(compile_check)) {
    list()
  } else {
    stop2("Unknown compile_check value")
  }
}

# fit Stan model
# @param model Stan model code
# @return validated Stan model code
fit_model <- function(model, backend, ...) {
  backend <- as_one_character(backend)
  .fit_model <- get(paste0(".fit_model_", backend), mode = "function")
  .fit_model(model, ...)
}

# fit Stan model with rstan
# @param model a compiled Stan model
# @param sdata named list to be passed to Stan as data
# @return a fitted Stan model
.fit_model_rstan <- function(model, sdata, algorithm, iter, warmup, thin, 
                             chains, cores, threads, inits, exclude, seed, 
                             control, silent, future, ...) {
  
  # some input checks and housekeeping
  if (use_threading(threads)) {
    if (utils::packageVersion("rstan") >= 2.26) {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(rstan::rstan_options(threads_per_chain = threads_per_chain_def))
      rstan::rstan_options(threads_per_chain = threads$threads)
    } else {
      stop2("Threading is not supported by backend 'rstan' version ",
            utils::packageVersion("rstan"), ".")
    }
  }
  if (is.character(inits) && !inits %in% c("random", "0")) {
    inits <- get(inits, mode = "function", envir = parent.frame())
  }
  args <- nlist(
    object = model, data = sdata, iter, seed, 
    init = inits, pars = exclude, include = FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  
  # do the actual sampling
  if (silent < 2) {
    message("Start sampling") 
  }
  if (algorithm %in% c("sampling", "fixed_param")) {
    c(args) <- nlist(warmup, thin, control, show_messages = !silent)
    if (algorithm == "fixed_param") {
      args$algorithm <- "Fixed_param"
    }
    if (future) {
      if (cores > 1L) {
        warning2("Argument 'cores' is ignored when using 'future'.")
      }
      args$chains <- 1L
      futures <- fits <- vector("list", chains)
      for (i in seq_len(chains)) {
        args$chain_id <- i
        if (is.list(inits)) {
          args$init <- inits[i]
        }
        futures[[i]] <- future::future(
          brms::do_call(rstan::sampling, args), 
          packages = "rstan",
          seed = TRUE
        )
      }
      for (i in seq_len(chains)) {
        fits[[i]] <- future::value(futures[[i]]) 
      }
      out <- rstan::sflist2stanfit(fits)
      rm(futures, fits)
    } else {
      c(args) <- nlist(chains, cores)
      out <- do_call(rstan::sampling, args) 
    }
  } else if (algorithm %in% c("fullrank", "meanfield")) {
    # vb does not support parallel execution
    c(args) <- nlist(algorithm)
    out <- do_call(rstan::vb, args)
  } else {
    stop2("Algorithm '", algorithm, "' is not supported.")
  }
  out
}

# fit Stan model with cmdstanr
# @param model a compiled Stan model
# @param sdata named list to be passed to Stan as data
# @return a fitted Stan model
.fit_model_cmdstanr <- function(model, sdata, algorithm, iter, warmup, thin, 
                                chains, cores, threads, inits, exclude, seed, 
                                control, silent, future, ...) {
  
  require_package("cmdstanr")
  # some input checks and housekeeping
  class(sdata) <- "list"
  if (isNA(seed)) {
    seed <- NULL
  }
  if (is_equal(inits, "random")) {
    inits <- NULL
  } else if (is_equal(inits, "0")) {
    inits <- 0
  }
  if (future) {
    stop2("Argument 'future' is not supported by backend 'cmdstanr'.")
  }
  args <- nlist(data = sdata, seed, init = inits)
  if (use_threading(threads)) {
    args$threads_per_chain <- threads$threads
  }
  # TODO: exclude variables via 'exclude'
  dots <- list(...)
  args[names(dots)] <- dots
  args[names(control)] <- control
  
  chains <- as_one_numeric(chains)
  empty_model <- chains <= 0
  if (empty_model) {
    # fit the model with minimal amount of draws
    # TODO: replace with a better solution
    chains <- 1
    iter <- 2
    warmup <- 1
    thin <- 1
    cores <- 1
  }
  
  # do the actual sampling
  if (silent < 2) {
    message("Start sampling") 
  }
  if (algorithm %in% c("sampling", "fixed_param")) {
    c(args) <- nlist(
      iter_sampling = iter - warmup,
      iter_warmup = warmup, 
      chains, thin, 
      parallel_chains = cores,
      show_messages = !silent,
      fixed_param = algorithm == "fixed_param"
    )
    out <- do_call(model$sample, args)
  } else if (algorithm %in% c("fullrank", "meanfield")) {
    # vb does not support parallel execution
    c(args) <- nlist(iter, algorithm)
    out <- do_call(model$variational, args)
  } else {
    stop2("Algorithm '", algorithm, "' is not supported.")
  }
  # transform into stanfit object for consistent output structure
  out <- rstan::read_stan_csv(out$output_files())
  # allow updating the model without recompilation
  attributes(out)$CmdStanModel <- model
  if (empty_model) {
    # allow correct updating of an 'empty' model
    out@sim <- list()
  }
  out
}


.fit_model_mock <- function(model, sdata, algorithm, iter, warmup, thin, 
chains, cores, threads, inits, exclude, seed, 
control, silent, future, stanfit, ...) {
  if(is.function(stanfit)) {
    stanfit()
  } else {
    stanfit
  }
}

# extract the compiled model
# @param x brmsfit object
compiled_model <- function(x) {
  stopifnot(is.brmsfit(x))
  backend <- x$backend %||% "rstan"
  if (backend == "rstan") {
    out <- rstan::get_stanmodel(x$fit)
  } else if (backend == "cmdstanr") {
    out <- attributes(x$fit)$CmdStanModel
  } else if(backend == "mock") {
    stop2("Compiled models not supported in mock backend")
  }
  out
}

# supported Stan backends
backend_choices <- function() {
  c("rstan", "cmdstanr", "mock")
}

# supported Stan algorithms
algorithm_choices <- function() {
  c("sampling", "meanfield", "fullrank", "fixed_param")
}

# check if the model was fit the the required backend
require_backend <- function(backend, x) {
  stopifnot(is.brmsfit(x))
  backend <- match.arg(backend, backend_choices())
  if (isTRUE(x$backend != backend)) {
    stop2("Backend '", backend, "' is required for this method.")
  }
  invisible(TRUE)
}

#' Threading in Stan
#' 
#' Use threads for within-chain parallelization in \pkg{Stan} via the \pkg{brms}
#' interface. Within-chain parallelization is experimental! We recommend its use
#' only if you are experienced with Stan's \code{reduce_sum} function and have a
#' slow running model that cannot be sped up by any other means.
#' 
#' @param threads Number of threads to use in within-chain parallelization.
#' @param grainsize Number of observations evaluated together in one chunk on
#'   one of the CPUs used for threading. If \code{NULL} (the default),
#'   \code{grainsize} is currently chosen as \code{max(100, N / (2 *
#'   threads))}, where \code{N} is the number of observations in the data. This
#'   default is experimental and may change in the future without prior notice.
#' @param static Logical. Apply the static (non-adaptive) version of
#'   \code{reduce_sum}? Defaults to \code{FALSE}. Setting it to \code{TRUE}
#'   is required to achieve exact reproducibility of the model results
#'   (if the random seed is set as well).
#' 
#' @return A \code{brmsthreads} object which can be passed to the
#'   \code{threads} argument of \code{brm} and related functions.
#'
#' @details The adaptive scheduling procedure used by \code{reduce_sum} will
#'   prevent the results to be exactly reproducible even if you set the random
#'   seed. If you need exact reproducibility, you have to set argument
#'   \code{static = TRUE} which may reduce efficiency a bit.
#'
#'   To ensure that chunks (whose size is defined by \code{grainsize}) require
#'   roughly the same amount of computing time, we recommend storing
#'   observations in random order in the data. At least, please avoid sorting
#'   observations after the response values. This is because the latter often
#'   cause variations in the computing time of the pointwise log-likelihood, 
#'   which makes up a big part of the parallelized code.
#'   
#' @examples 
#' \dontrun{
#' # this model just serves as an illustration
#' # threading may not actually speed things up here
#' fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'            data = epilepsy, family = negbinomial(),
#'            chains = 1, threads = threading(2, grainsize = 100),
#'            backend = "cmdstanr")
#' summary(fit)
#' }
#' 
#' @export
threading <- function(threads = NULL, grainsize = NULL, static = FALSE) {
  out <- list(threads = NULL, grainsize = NULL)
  class(out) <- "brmsthreads"
  if (!is.null(threads)) {
    threads <- as_one_numeric(threads)
    if (!is_wholenumber(threads) || threads < 1) {
      stop2("Number of threads needs to be positive.")
    }
    out$threads <- threads
  }
  if (!is.null(grainsize)) {
    grainsize <- as_one_numeric(grainsize)
    if (!is_wholenumber(grainsize) || grainsize < 1) {
      stop2("The grainsize needs to be positive.")
    }
    out$grainsize <- grainsize
  }
  out$static <- as_one_logical(static)
  out
}

is.brmsthreads <- function(x) {
  inherits(x, "brmsthreads")
}

# validate 'thread' argument
validate_threads <- function(threads) {
  if (is.null(threads)) {
    threads <- threading()
  } else if (is.numeric(threads)) {
    threads <- as_one_numeric(threads)
    threads <- threading(threads)
  } else if (!is.brmsthreads(threads)) {
    stop2("Argument 'threads' needs to be numeric or ", 
          "specified via the 'threading' function.")
  }
  threads
}

# is threading activated?
use_threading <- function(threads) {
  isTRUE(validate_threads(threads)$threads > 0)
}

# validate the 'silent' argument
validate_silent <- function(silent) {
  silent <- as_one_integer(silent)
  if (silent < 0 || silent > 2) {
    stop2("'silent' must be between 0 and 2.")
  }
  silent
}
