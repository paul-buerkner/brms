# Helper function that will capture formal elements of brm function and return
# as a list
brm_to_call <- function(formula, data = NULL, family = gaussian(), prior = NULL,
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
                        file_auto = getOption("brms.file_auto", FALSE),
                        empty = FALSE, rename = TRUE, ...) {

  # validate arguments later passed to Stan
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
  file_auto<- as_one_logical(file_auto)

  .call <- match.call()
  orig_seed <- seed
  # This list must include only/all the parameters that may change the result
  args_list <- nlist(formula, data, family, prior, autocor, data2, cov_ranef,
                     sample_prior, sparse, knots, drop_unused_levels, stanvars,
                     stan_funs, fit, save_pars, save_ranef, save_mevars,
                     save_all_pars, init, inits, chains, iter, warmup, thin,
                     cores, threads, opencl, normalize, control, algorithm,
                     backend, future, orig_seed, stan_model_args, empty , .call)
  args_list
}



## Recursively sanitize brms inputs for hashing purposes
## Converts or strips elements that may cause unnecessary variability in hashing
sanitizer_recursive <- function(x) {

  # Handle formula objects by removing their environments and converting to character
  if (inherits(x, "formula")) {
    environment(x) <- emptyenv()
    return(as.character(x))
  }

  # Handle brmsformula objects by recursively sanitizing internal formulas
  if (inherits(x, "brmsformula")) {
    x$formula <- sanitizer_recursive(x$formula)
    if (!is.null(x$pforms)) x$pforms <- lapply(x$pforms, sanitizer_recursive)
    if (!is.null(x$nlpars)) x$nlpars <- sort(x$nlpars)  # Sort for consistent ordering
    return(x)
  }

  # Handle multivariate brms formulas (mvbrmsformula) by sorting and sanitizing each component
  if (inherits(x, "mvbrmsformula")) {
    # Sort component names to ensure ordering does not affect the hash
    x$forms <- lapply(x$forms[sort(names(x$forms))], sanitizer_recursive)
    return(x)
  }

  # Handle family objects by reducing to essential components
  if (inherits(x, "family")) {
    return(list(family = x$family, link = x$link))
  }

  # Handle function objects by hashing their body as character
  if (is.function(x)) {
    return(deparse(body(x), width.cutoff = 500L))
  }

  # Handle calls, language objects, and expressions by converting to character
  if (is.call(x) || is.language(x) || is.expression(x)) {
    return(deparse(x, width.cutoff = 500L))
  }

  # Recursively handle generic lists (excluding data frames), sorting for stable order
  if (is.list(x) && !inherits(x, "data.frame")) {
    x <- x[sort(names(x))]
    return(lapply(x, sanitizer_recursive))
  }

  # Return all other values unchanged
  x
}

# Clean and sanitize the formula object for hashing
handle_formula_for_hash<-function(formula){
  sanitizer_recursive(formula)
}
# Clean and sanitize the family object for hashing
handle_family_for_hash<-function(family){
  sanitizer_recursive(family)
}
# Generate a unique hash string based on key components of the brm() call
hash_brm_call_list <-function(args_list){
  algo = "xxhash64"
  requireNamespace("digest")
  data <- args_list$data
  # handle formula
  formula <- handle_formula_for_hash(args_list$formula)
  # handle family
  family <- handle_family_for_hash(args_list$family)
  # handle data
  data <- get_same_data_if_null(data)
  hashed_data  <- hash_data(data)
  #package version
  brms_version <- packageVersion("brms")
  args_list <- nlist(formula, family, hashed_data, brms_version )
  digest::digest(args_list, algo = algo)
}

# Internal function to use testing that will behave like brm function
# collect the arguments but return hash value for the call
hash_func_for_testing <- function(...){

  # 1 - get list from brm()
  dots <- brm_to_call(...)
  # 2 - convert list to hash value
  hash <-  hash_brm_call_list(dots)
  # create file argument value
  hash

}

# If the file_auto argument is TRUE, generate a file name based on the model inputs
# to automatically save and reuse fitted model results.
# If file_auto is FALSE, return the original file and file_refit values unchanged.
create_filename_auto <- function(file, file_refit, file_auto, args_list) {
  if (!file_auto) {
    return(nlist(file, file_refit))
  }
  hash <- hash_brm_call_list(args_list)
  orig_file <- file
  file <- paste0('cache-brm-result_', hash, '.Rds')
  orig_file_refit <- file_refit
  file_refit <- "on_change"
  # We inform user that we override file or file_refit arguments in case necessary
  if (!is.null(orig_file) | orig_file_refit != 'on_change') {
    message("Since file_auto = TRUE, the file and file_refit arguments were overwritten.")
  }
  nlist(file, file_refit)
}

..f_testing <- function(){



  # same all
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
             data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~  zAge + zBase * Trt + (1|patient),
                     data = epilepsy, family = poisson() )

  expect_true( c1 == c2 )

  # different formula
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~   zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson() )

  expect_false( c1 == c2 )

  # different data
  c1 <-   hash_func_for_testing(count ~  zAge + zBase * Trt + (1|patient),
                                data = epilepsy, family = poisson())

  c2 <-  hash_func_for_testing( count ~   zBase * Trt + (1|patient),
                                data = epilepsy[-c(1), ], family = poisson() )

  expect_false( c1 == c2 )



  a <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  b <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  expect_true( a == b )

  a <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(2) , ], family = poisson() )

  b <- hash_func_for_testing( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  expect_false( a == b )



  m1 <- brm( count ~  zBase * Trt + (1|patient),
             data = epilepsy[-c(2) , ], family = poisson()  , file = "m")


}





