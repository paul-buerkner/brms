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
.parse_model_rstan <- function(model, silent = TRUE, ...) {
  out <- eval_silent(
    rstan::stanc(model_code = model, ...),
    type = "message", 
    try = TRUE, 
    silent = silent
  )
  out$model_code
}

# parse Stan model code with cmdstanr
# @param model Stan model code
# @return validated Stan model code
.parse_model_cmdstanr <- function(model, silent = TRUE, ...) {
  require_package("cmdstanr")
  temp_file <- cmdstanr::write_stan_tempfile(model)
  out <- eval_silent(
    cmdstanr::cmdstan_model(temp_file, compile = FALSE, ...),
    type = "message", 
    try = TRUE, 
    silent = silent
  )
  collapse(out$code(), "\n")
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
.compile_model_rstan <- function(model, ...) {
  args <- list(...)
  silent <- !length(args)
  args$model_code <- model
  message("Compiling Stan program...")
  eval_silent(
    do_call(rstan::stan_model, args),
    silent = silent, type = "message"
  )
}

# compile Stan model with cmdstanr
# @param model Stan model code
# @return model compiled with cmdstanr
.compile_model_cmdstanr <- function(model, ...) {
  require_package("cmdstanr")
  args <- list(...)
  silent <- !length(args)
  args$stan_file <- cmdstanr::write_stan_tempfile(model)
  eval_silent(
    do_call(cmdstanr::cmdstan_model, args),
    silent = silent, type = "message"
  )
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
                             chains, cores, inits, exclude, seed, control, 
                             silent, future, ...) {
  
  # some input checks and housekeeping
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
  message("Start sampling")
  if (algorithm == "sampling") {
    c(args) <- nlist(warmup, thin, control, show_messages = !silent)
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
          packages = "rstan"
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
                                chains, cores, inits, exclude, seed, control, 
                                silent, future, ...) {
  
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
  # TODO: exclude variables via 'exclude'
  # TODO: silence messages via 'silent'
  dots <- list(...)
  args[names(dots)] <- dots
  args[names(control)] <- control
  
  # do the actual sampling
  message("Start sampling")
  if (algorithm == "sampling") {
    c(args) <- nlist(
      iter_sampling = iter - warmup,
      iter_warmup = warmup, 
      chains, cores, thin
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
  out
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
  }
  out
}

# supported Stan backends
backend_choices <- function() {
  c("rstan", "cmdstanr")
}

# supported Stan algorithms
algorithm_choices <- function() {
  c("sampling", "meanfield", "fullrank")
}
