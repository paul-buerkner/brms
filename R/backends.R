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
  # if (cmdstanr::cmdstan_version() >= "2.29.0") {
  #   .canonicalize_stan_model(temp_file, overwrite_file = TRUE)
  # }
  out <- eval_silent(
    cmdstanr::cmdstan_model(temp_file, compile = FALSE, ...),
    type = "message", try = TRUE, silent = silent
  )
  out$check_syntax(quiet = TRUE)
  collapse(out$code(), "\n")
}

# parse model with a mock backend for testing
.parse_model_mock <- function(model, silent = TRUE, parse_error = NULL,
                              parse_check = "rstan", ...) {
  if (!is.null(parse_error)) {
    stop2(parse_error)
  } else if (parse_check == "rstan") {
    out <- .parse_model_rstan(model, silent = silent, ...)
  } else if (parse_check == "cmdstanr") {
    out <- .parse_model_cmdstanr(model, silent = silent, ...)
  } else if (is.null(parse_check)) {
    out <- "mock_code"
  } else {
    stop2("Unknown 'parse_check' value.")
  }
  out
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
.compile_model_rstan <- function(model, threads, opencl, silent = 1, ...) {
  args <- list(...)
  args$model_code <- model
  if (silent < 2) {
    message("Compiling Stan program...")
  }
  if (use_threading(threads)) {
    if (utils::packageVersion("rstan") >= "2.26") {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(rstan::rstan_options(threads_per_chain = threads_per_chain_def))
      rstan::rstan_options(threads_per_chain = threads$threads)
    } else {
      stop2("Threading is not supported by backend 'rstan' version ",
            utils::packageVersion("rstan"), ".")
    }
  }
  if (use_opencl(opencl)) {
    stop2("OpenCL is not supported by backend 'rstan' version ",
          utils::packageVersion("rstan"), ".")
  }
  eval_silent(
    do_call(rstan::stan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

# compile Stan model with cmdstanr
# @param model Stan model code
# @return model compiled with cmdstanr
.compile_model_cmdstanr <- function(model, threads, opencl, silent = 1, ...) {
  require_package("cmdstanr")
  args <- list(...)
  args$stan_file <- cmdstanr::write_stan_file(model)
  # if (cmdstanr::cmdstan_version() >= "2.29.0") {
  #   .canonicalize_stan_model(args$stan_file, overwrite_file = TRUE)
  # }
  if (use_threading(threads)) {
    args$cpp_options$stan_threads <- TRUE
  }
  if (use_opencl(opencl)) {
    args$cpp_options$stan_opencl <- TRUE
  }
  eval_silent(
    do_call(cmdstanr::cmdstan_model, args),
    type = "message", try = TRUE, silent = silent >= 2
  )
}

# compile model with a mock backend for testing
.compile_model_mock <- function(model, threads, opencl, compile_check = "rstan",
                                compile_error = NULL, silent = 1, ...) {
  if (!is.null(compile_error)) {
    stop2(compile_error)
  } else if (compile_check == "rstan") {
    out <- .parse_model_rstan(model, silent = silent, ...)
  } else if (compile_check == "cmdstanr") {
    out <- .parse_model_cmdstanr(model, silent = silent, ...)
  } else if (is.null(compile_check)) {
    out <- list()
  } else {
    stop2("Unknown 'compile_check' value.")
  }
  out
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
                             chains, cores, threads, opencl, init, exclude,
                             seed, control, silent, future, ...) {

  # some input checks and housekeeping
  if (use_threading(threads)) {
    if (utils::packageVersion("rstan") >= "2.26") {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(rstan::rstan_options(threads_per_chain = threads_per_chain_def))
      rstan::rstan_options(threads_per_chain = threads$threads)
    } else {
      stop2("Threading is not supported by backend 'rstan' version ",
            utils::packageVersion("rstan"), ".")
    }
  }
  if (use_opencl(opencl)) {
    stop2("OpenCL is not supported by backend 'rstan' version ",
          utils::packageVersion("rstan"), ".")
  }
  if (is.null(init)) {
    init <- "random"
  } else if (is.character(init) && !init %in% c("random", "0")) {
    init <- get(init, mode = "function", envir = parent.frame())
  }
  args <- nlist(
    object = model, data = sdata, iter, seed,
    init = init, pars = exclude, include = FALSE
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
        if (is.list(init)) {
          args$init <- init[i]
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
  out <- repair_stanfit(out)
  out
}

# fit Stan model with cmdstanr
# @param model a compiled Stan model
# @param sdata named list to be passed to Stan as data
# @return a fitted Stan model
.fit_model_cmdstanr <- function(model, sdata, algorithm, iter, warmup, thin,
                                chains, cores, threads, opencl, init, exclude,
                                seed, control, silent, future, ...) {

  require_package("cmdstanr")
  # some input checks and housekeeping
  class(sdata) <- "list"
  if (isNA(seed)) {
    seed <- NULL
  }
  if (is_equal(init, "random")) {
    init <- NULL
  } else if (is_equal(init, "0")) {
    init <- 0
  }
  if (future) {
    stop2("Argument 'future' is not supported by backend 'cmdstanr'.")
  }
  args <- nlist(data = sdata, seed, init)
  if (use_threading(threads)) {
    if (algorithm %in% c("sampling", "fixed_param")) {
      args$threads_per_chain <- threads$threads
    } else if (algorithm %in% c("fullrank", "meanfield")) {
      args$threads <- threads$threads
    }
  }
  if (use_opencl(opencl)) {
    args$opencl_ids <- opencl$ids
  }
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
      show_messages = silent < 2,
      show_exceptions = silent == 0,
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

  out <- read_csv_as_stanfit(
    out$output_files(), variables = out$metadata()$variables,
    model = model, exclude = exclude
  )

  if (empty_model) {
    # allow correct updating of an 'empty' model
    out@sim <- list()
  }
  out
}

# fit model with a mock backend for testing
.fit_model_mock <- function(model, sdata, algorithm, iter, warmup, thin,
                            chains, cores, threads, opencl, init, exclude,
                            seed, control, silent, future, mock_fit, ...) {
  if (is.function(mock_fit)) {
    out <- mock_fit()
  } else {
    out <- mock_fit
  }
  out
}

# extract the compiled stan model
# @param x brmsfit object
compiled_model <- function(x) {
  stopifnot(is.brmsfit(x))
  backend <- x$backend %||% "rstan"
  if (backend == "rstan") {
    out <- rstan::get_stanmodel(x$fit)
  } else if (backend == "cmdstanr") {
    out <- attributes(x$fit)$CmdStanModel
  } else if (backend == "mock") {
    stop2("'compiled_model' is not supported in the mock backend.")
  }
  out
}

# Does the model need recompilation before being able to sample again?
needs_recompilation <- function(x) {
  stopifnot(is.brmsfit(x))
  backend <- x$backend %||% "rstan"
  if (backend == "rstan") {
    # TODO: figure out when rstan requires recompilation
    out <- FALSE
  } else if (backend == "cmdstanr") {
    exe_file <- attributes(x$fit)$CmdStanModel$exe_file()
    out <- !is.character(exe_file) || !file.exists(exe_file)
  } else if (backend == "mock") {
    out <- FALSE
  }
  out
}

#' Recompile Stan models in \code{brmsfit} objects
#'
#' Recompile the Stan model inside a \code{brmsfit} object, if necessary.
#' This does not change the model, it simply recreates the executable
#' so that sampling is possible again.
#'
#' @param x An object of class \code{brmsfit}.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled. If \code{NULL} (the default), \code{recompile_model} tries
#'   to figure out internally, if recompilation is necessary. Setting it to
#'   \code{FALSE} will cause \code{recompile_model} to always return the
#'   \code{brmsfit} object unchanged.
#'
#' @return A (possibly updated) \code{brmsfit} object.
#'
#' @export
recompile_model <- function(x, recompile = NULL) {
  stopifnot(is.brmsfit(x))
  if (is.null(recompile)) {
    recompile <- needs_recompilation(x)
  }
  recompile <- as_one_logical(recompile)
  if (!recompile) {
    return(x)
  }
  message("Recompiling the Stan model")
  backend <- x$backend %||% "rstan"
  new_model <- compile_model(
    stancode(x), backend = backend, threads = x$threads,
    opencl = x$opencl, silent = 2
  )
  if (backend == "rstan") {
    x$fit@stanmodel <- new_model
  } else if (backend == "cmdstanr") {
    attributes(x)$CmdStanModel <- new_model
  } else if (backend == "mock") {
    stop2("'recompile_model' is not supported in the mock backend.")
  }
  x
}

# extract the elapsed time during model fitting
# @param x brmsfit object
elapsed_time <- function(x) {
  stopifnot(is.brmsfit(x))
  backend <- x$backend %||% "rstan"
  if (backend == "rstan") {
    out <- rstan::get_elapsed_time(x$fit)
    out <- data.frame(
      chain_id = seq_len(nrow(out)),
      warmup = out[, "warmup"],
      sampling = out[, "sample"]
    )
    out$total <- out$warmup + out$sampling
    rownames(out) <- NULL
  } else if (backend == "cmdstanr") {
    out <- attributes(x$fit)$metadata$time$chains
  } else if (backend == "mock") {
    stop2("'elapsed_time' not supported in the mock backend.")
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

#' GPU support in Stan via OpenCL
#'
#' Use OpenCL for GPU support in \pkg{Stan} via the \pkg{brms} interface. Only
#' some \pkg{Stan} functions can be run on a GPU at this point and so
#' a lot of \pkg{brms} models won't benefit from OpenCL for now.
#'
#' @param ids (integer vector of length 2) The platform and device IDs of the
#'   OpenCL device to use for fitting. If you don't know the IDs of your OpenCL
#'   device, \code{c(0,0)} is most likely what you need.
#'
#' @return A \code{brmsopencl} object which can be passed to the
#'   \code{opencl} argument of \code{brm} and related functions.
#'
#' @details For more details on OpenCL in \pkg{Stan}, check out
#' \url{https://mc-stan.org/docs/2_26/cmdstan-guide/parallelization.html#opencl}
#' as well as \url{https://mc-stan.org/docs/2_26/stan-users-guide/opencl.html}.
#'
#' @examples
#' \dontrun{
#' # this model just serves as an illustration
#' # OpenCL may not actually speed things up here
#' fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'            data = epilepsy, family = poisson(),
#'            chains = 2, cores = 2, opencl = opencl(c(0, 0)),
#'            backend = "cmdstanr")
#' summary(fit)
#' }
#'
#' @export
opencl <- function(ids = NULL) {
  out <- list(ids = NULL)
  class(out) <- "brmsopencl"
  if (!is.null(ids)) {
    ids <- as.integer(ids)
    if (!length(ids) == 2L) {
      stop2("OpenCl 'ids' needs to be an integer vector of length 2.")
    }
    out$ids <- ids
  }
  out
}

is.brmsopencl <- function(x) {
  inherits(x, "brmsopencl")
}

# validate the 'opencl' argument
validate_opencl <- function(opencl) {
  if (is.null(opencl)) {
    opencl <- opencl()
  } else if (is.numeric(opencl)) {
    opencl <- opencl(opencl)
  } else if (!is.brmsopencl(opencl)) {
    stop2("Argument 'opencl' needs to an integer vector or ",
          "specified via the 'opencl' function.")
  }
  opencl
}

# is OpenCL activated?
use_opencl <- function(opencl) {
  !is.null(validate_opencl(opencl)$ids)
}

# validate the 'silent' argument
validate_silent <- function(silent) {
  silent <- as_one_integer(silent)
  if (silent < 0 || silent > 2) {
    stop2("'silent' must be between 0 and 2.")
  }
  silent
}

# ensure that variable dimensions at the end are correctly written
# convert names like b.1.1 to b[1,1]
repair_variable_names <- function(x) {
  x <- sub("\\.", "[", x)
  x <- gsub("\\.", ",", x)
  x[grep("\\[", x)] <- paste0(x[grep("\\[", x)], "]")
  x
}

# repair parameter names of stanfit objects
repair_stanfit <- function(x) {
  stopifnot(is.stanfit(x))
  if (!length(x@sim$fnames_oi)) {
    # nothing to rename
    return(x)
  }
  # the posterior package cannot deal with non-unique parameter names
  # this case happens rarely but might happen when sample_prior = "yes"
  x@sim$fnames_oi <- make.unique(as.character(x@sim$fnames_oi), "__")
  for (i in seq_along(x@sim$samples)) {
    # stanfit may have renamed dimension suffixes (#1218)
    if (length(x@sim$samples[[i]]) == length(x@sim$fnames_oi)) {
      names(x@sim$samples[[i]]) <- x@sim$fnames_oi
    }
  }
  x
}

# possible options for argument 'file_refit'
file_refit_options <- function() {
  c("never", "always", "on_change")
}

# canonicalize Stan model file in accordance with the current Stan version
# this function may no longer be needed due to rstan 2.26+ now being on CRAN
# for more details see https://github.com/paul-buerkner/brms/issues/1544
# .canonicalize_stan_model <- function(stan_file, overwrite_file = TRUE) {
#   cmdstan_mod <- cmdstanr::cmdstan_model(stan_file, compile = FALSE)
#   out <- utils::capture.output(
#     cmdstan_mod$format(
#       canonicalize = list("deprecations", "braces", "parentheses"),
#       overwrite_file = overwrite_file, backup = FALSE
#     )
#   )
#   paste0(out, collapse = "\n")
# }

#' Read CmdStan CSV files as a brms-formatted stanfit object
#'
#' \code{read_csv_as_stanfit} is used internally to read CmdStan CSV files into a
#' \code{stanfit} object that is consistent with the structure of the fit slot of a
#' brmsfit object.
#'
#' @param files Character vector of CSV files names where draws are stored.
#' @param variables Character vector of variables to extract from the CSV files.
#' @param sampler_diagnostics Character vector of sampler diagnostics to extract.
#' @param model A compiled cmdstanr model object (optional). Provide this argument
#'  if you want to allow updating the model without recompilation.
#' @param exclude Character vector of variables to exclude from the stanfit. Only
#'  used when \code{variables} is also specified.
#'
#' @return A stanfit object consistent with the structure of the \code{fit}
#'  slot of a brmsfit object.
#'
#' @examples
#' \dontrun{
#' # fit a model manually via cmdstanr
#' scode <- stancode(count ~ Trt, data = epilepsy)
#' sdata <- standata(count ~ Trt, data = epilepsy)
#' mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(scode))
#' stanfit <- mod$sample(data = sdata)
#'
#' # feed the Stan model back into brms
#' fit <- brm(count ~ Trt, data = epilepsy, empty = TRUE, backend = 'cmdstanr')
#' fit$fit <- read_csv_as_stanfit(stanfit$output_files(), model = mod)
#' fit <- rename_pars(fit)
#' summary(fit)
#' }
#'
#' @export
read_csv_as_stanfit <- function(files, variables = NULL, sampler_diagnostics = NULL,
                                model = NULL, exclude = "") {
  require_package("cmdstanr")

  if (!is.null(variables)) {
    # ensure that only relevant variables are read from CSV
    variables <- repair_variable_names(variables)
    variables <- unique(sub("\\[.+", "", variables))
    variables <- setdiff(variables, exclude)
    # temp fix for cmdstanr not recognizing the variable names it produces  #1473
    variables <- ifelse(variables == "lp_approx__", "log_g__", variables)
  }

  csfit <- cmdstanr::read_cmdstan_csv(
    files = files, variables = variables,
    sampler_diagnostics = sampler_diagnostics,
    format = NULL
  )

  # @model_name
  model_name = gsub(".csv", "", basename(files[[1]]))

  # @model_pars
  svars <- variables %||% csfit$metadata$stan_variables
  if ("lp__" %in% svars) {
    svars <- c(setdiff(svars, "lp__"), "lp__")
  }
  pars_oi <- svars
  par_names <- csfit$metadata$model_params

  # @par_dims
  par_dims <- vector("list", length(svars))

  names(par_dims) <- svars
  par_dims <- lapply(par_dims, function(x) x <- integer(0))

  pdims_num <- ulapply(
    svars, function(x) sum(grepl(paste0("^", x, "\\[.*\\]$"), par_names))
  )
  par_dims[pdims_num != 0] <-
    csfit$metadata$stan_variable_sizes[svars][pdims_num != 0]

  # @mode
  mode <- 0L

  # @sim
  rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                         "divergent__", "n_leapfrog__", "energy__")

  if (!is.null(sampler_diagnostics)) {
    rstan_diagn_order <- rstan_diagn_order[rstan_diagn_order %in% sampler_diagnostics]
  }

  res_vars <- c(".chain", ".iteration", ".draw")
  if ("post_warmup_draws" %in% names(csfit)) {
    # for MCMC samplers
    n_chains <- max(
      nchains(csfit$warmup_draws),
      nchains(csfit$post_warmup_draws)
    )
    n_iter_warmup <- niterations(csfit$warmup_draws)
    n_iter_sample <- niterations(csfit$post_warmup_draws)
    if (n_iter_warmup > 0) {
      csfit$warmup_draws <- as_draws_df(csfit$warmup_draws)
      csfit$warmup_sampler_diagnostics <-
        as_draws_df(csfit$warmup_sampler_diagnostics)
    }
    if (n_iter_sample > 0) {
      csfit$post_warmup_draws <- as_draws_df(csfit$post_warmup_draws)
      csfit$post_warmup_sampler_diagnostics <-
        as_draws_df(csfit$post_warmup_sampler_diagnostics)
    }

    # called 'samples' for consistency with rstan
    samples <- rbind(csfit$warmup_draws, csfit$post_warmup_draws)
    # manage memory
    csfit$warmup_draws <- NULL
    csfit$post_warmup_draws <- NULL

    # prepare sampler diagnostics
    diagnostics <- rbind(csfit$warmup_sampler_diagnostics,
                         csfit$post_warmup_sampler_diagnostics)
    # manage memory
    csfit$warmup_sampler_diagnostics <- NULL
    csfit$post_warmup_sampler_diagnostics <- NULL
    # convert to regular data.frame
    diagnostics <- as.data.frame(diagnostics)
    diag_chain_ids <- diagnostics$.chain
    diagnostics[res_vars] <- NULL

  } else if ("draws" %in% names(csfit)) {
    # for variational inference "samplers"
    n_chains <- 1
    n_iter_warmup <- 0
    n_iter_sample <- niterations(csfit$draws)
    if (n_iter_sample > 0) {
      csfit$draws <- as_draws_df(csfit$draws)
    }

    # called 'samples' for consistency with rstan
    samples <- csfit$draws
    # manage memory
    csfit$draws <- NULL

    # VI has no sampler diagnostics
    diag_chain_ids <- rep(1L, nrow(samples))
    diagnostics <- as.data.frame(matrix(nrow = nrow(samples), ncol = 0))
  }

  # convert to regular data.frame
  samples <- as.data.frame(samples)
  chain_ids <- samples$.chain
  samples[res_vars] <- NULL
  if ("lp__" %in% colnames(samples)) {
    samples <- move2end(samples, "lp__")
  }

  fnames_oi <- colnames(samples)

  # split samples into chains
  samples <- split(samples, chain_ids)
  names(samples) <- NULL

  # split diagnostics into chains
  diagnostics <- split(diagnostics, diag_chain_ids)
  names(diagnostics) <- NULL

  #  @sim$sample: largely 113-130 from rstan::read_stan_csv
  values <- list()
  values$algorithm <- csfit$metadata$algorithm
  values$engine <- csfit$metadata$engine
  values$metric <- csfit$metadata$metric

  sampler_t <- NULL
  if (!is.null(values$algorithm)) {
    if (values$algorithm == "rwm" || values$algorithm == "Metropolis") {
      sampler_t <- "Metropolis"
    } else if (values$algorithm == "hmc") {
      if (values$engine == "static") {
        sampler_t <- "HMC"
      } else {
        if (values$metric == "unit_e") {
          sampler_t <- "NUTS(unit_e)"
        } else if (values$metric == "diag_e") {
          sampler_t <- "NUTS(diag_e)"
        } else if (values$metric == "dense_e") {
          sampler_t <- "NUTS(dense_e)"
        }
      }
    }
  }

  adapt_info <- vector("list", 4)
  idx_samples <- (n_iter_warmup + 1):(n_iter_warmup + n_iter_sample)

  for (i in seq_along(samples)) {
    m <- colMeans(samples[[i]][idx_samples, , drop = FALSE])
    rownames(samples[[i]]) <- seq_rows(samples[[i]])
    attr(samples[[i]], "sampler_params") <- diagnostics[[i]][rstan_diagn_order]
    rownames(attr(samples[[i]], "sampler_params")) <- seq_rows(diagnostics[[i]])

    # reformat back to text
    if (length(csfit$inv_metric)) {
      if (is_equal(sampler_t, "NUTS(dense_e)")) {
        mmatrix_txt <- "\n# Elements of inverse mass matrix:\n# "
        mmat <- paste0(apply(csfit$inv_metric[[i]], 1, paste0, collapse = ", "),
                       collapse = "\n# ")
      } else {
        mmatrix_txt <- "\n# Diagonal elements of inverse mass matrix:\n# "
        mmat <- paste0(csfit$inv_metric[[i]], collapse = ", ")
      }

      adapt_info[[i]] <- paste0("# Step size = ",
                                csfit$step_size[[i]],
                                mmatrix_txt,
                                mmat, "\n# ")
      attr(samples[[i]], "adaptation_info") <- adapt_info[[i]]
    } else {
      attr(samples[[i]], "adaptation_info") <- character(0)
    }


    attr(samples[[i]], "args") <- list(sampler_t = sampler_t, chain_id = i)

    if (NROW(csfit$metadata$time)) {
      time_i <- as.double(csfit$metadata$time[i, c("warmup", "sampling")])
      names(time_i) <- c("warmup", "sample")
      attr(samples[[i]], "elapsed_time") <- time_i
    }

    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }

  perm_lst <- lapply(seq_len(n_chains), function(id) sample.int(n_iter_sample))

  # @sim
  sim <- list(
    samples = samples,
    iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
    thin = csfit$metadata$thin,
    warmup = csfit$metadata$iter_warmup,
    chains = n_chains,
    n_save = rep(n_iter_sample + n_iter_warmup, n_chains),
    warmup2 = rep(n_iter_warmup, n_chains),
    permutation = perm_lst,
    pars_oi = pars_oi,
    dims_oi = par_dims,
    fnames_oi = fnames_oi,
    n_flatnames = length(fnames_oi)
  )

  # @stan_args
  sargs <- list(
    stan_version_major = as.character(csfit$metadata$stan_version_major),
    stan_version_minor = as.character(csfit$metadata$stan_version_minor),
    stan_version_patch = as.character(csfit$metadata$stan_version_patch),
    model = csfit$metadata$model_name,
    start_datetime = gsub(" ", "", csfit$metadata$start_datetime),
    method = csfit$metadata$method,
    iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
    warmup = csfit$metadata$iter_warmup,
    save_warmup = csfit$metadata$save_warmup,
    thin = csfit$metadata$thin,
    engaged = as.character(csfit$metadata$adapt_engaged),
    gamma = csfit$metadata$gamma,
    delta = csfit$metadata$adapt_delta,
    kappa = csfit$metadata$kappa,
    t0 = csfit$metadata$t0,
    init_buffer = as.character(csfit$metadata$init_buffer),
    term_buffer = as.character(csfit$metadata$term_buffer),
    window = as.character(csfit$metadata$window),
    algorithm = csfit$metadata$algorithm,
    engine = csfit$metadata$engine,
    max_depth = csfit$metadata$max_treedepth,
    metric = csfit$metadata$metric,
    metric_file = character(0), # not stored in metadata
    stepsize = NA, # add in loop
    stepsize_jitter = csfit$metadata$stepsize_jitter,
    num_chains = as.character(csfit$metadata$num_chains),
    chain_id = NA, # add in loop
    file = character(0), # not stored in metadata
    init = NA, # add in loop
    seed = as.character(csfit$metadata$seed),
    file = NA, # add in loop
    diagnostic_file = character(0), # not stored in metadata
    refresh = as.character(csfit$metadata$refresh),
    sig_figs = as.character(csfit$metadata$sig_figs),
    profile_file = csfit$metadata$profile_file,
    num_threads = as.character(csfit$metadata$threads_per_chain),
    stanc_version = gsub(" ", "", csfit$metadata$stanc_version),
    stancflags = character(0), # not stored in metadata
    adaptation_info = NA, # add in loop
    has_time = is.numeric(csfit$metadata$time$total),
    time_info = NA, # add in loop
    sampler_t = sampler_t
  )

  sargs_rep <- replicate(n_chains, sargs, simplify = FALSE)

  for (i in seq_along(sargs_rep)) {
    sargs_rep[[i]]$chain_id <- i
    sargs_rep[[i]]$stepsize <- csfit$metadata$step_size[i]
    sargs_rep[[i]]$init <- as.character(csfit$metadata$init[i])
    # two 'file' elements: select the second
    file_idx <- which(names(sargs_rep[[i]]) == "file")
    sargs_rep[[i]][[file_idx[2]]] <- files[[i]]

    sargs_rep[[i]]$adaptation_info <- adapt_info[[i]]

    if (NROW(csfit$metadata$time)) {
      sargs_rep[[i]]$time_info <- paste0(
        c("#  Elapsed Time: ", "#                ", "#                ", "# "),
        c(csfit$metadata$time[i, c("warmup", "sampling", "total")], ""),
        c(" seconds (Warm-up)", " seconds (Sampling)", " seconds (Total)", "")
      )
    }
  }

  # @stanmodel
  null_dso <- new(
    "cxxdso", sig = list(character(0)), dso_saved = FALSE,
    dso_filename = character(0), modulename = character(0),
    system = R.version$system, cxxflags = character(0),
    .CXXDSOMISC = new.env(parent = emptyenv())
  )
  null_sm <- new(
    "stanmodel", model_name = model_name, model_code = character(0),
    model_cpp = list(), dso = null_dso
  )

  # @date
  sdate <- do.call(max, lapply(files, function(csv) file.info(csv)$mtime))
  sdate <- format(sdate, "%a %b %d %X %Y")

  out <- new(
    "stanfit",
    model_name = model_name,
    model_pars = svars,
    par_dims = par_dims,
    mode = mode,
    sim = sim,
    inits = list(),
    stan_args = sargs_rep,
    stanmodel = null_sm,
    date = sdate,  # not the time of sampling
    .MISC = new.env(parent = emptyenv())
  )

  attributes(out)$metadata <- csfit
  attributes(out)$CmdStanModel <- model
  out
}
