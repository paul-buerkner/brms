
hash_brms_callOLD <- function(args_list,
                              algo        = "xxhash64",
                              data_policy = c("names", "signature", "hash", "none")) {

  require_package("digest")

  ## ── 0. *NEW*  Make top-level arg order irrelevant ──────────────────────
  args_list <- args_list[sort(names(args_list))]

  ## ── 1. Sanitiser (recursive) ───────────────────────────────────────────
  clean <- function(x) {
    if (inherits(x, "formula")) {
      environment(x) <- emptyenv()
      return(as.character(x))
    }

    if (inherits(x, "brmsformula")) {
      x$formula <- clean(x$formula)
      if (!is.null(x$pforms)) x$pforms <- lapply(x$pforms, clean)
      if (!is.null(x$nlpars)) x$nlpars <- sort(x$nlpars)
      return(x)
    }

    if (inherits(x, "mvbrmsformula")) {
      ## *NEW*  sort component names so mvbf() order is irrelevant
      x$forms <- lapply(x$forms[sort(names(x$forms))], clean)
      return(x)
    }

    if (inherits(x, "family")) {
      return(list(family = x$family, link = x$link))
    }

    if (is.function(x)) {
      return(deparse(body(x), width.cutoff = 500L))
    }

    if (is.call(x) || is.language(x) || is.expression(x)) {
      return(deparse(x, width.cutoff = 500L))
    }

    if (is.list(x) && !inherits(x, "data.frame")) {
      x <- x[sort(names(x))]
      return(lapply(x, clean))
    }

    x
  }

  ## 2–5 unchanged ---------------------------------------------------------
  args_list <- lapply(args_list, clean)

  policy <- match.arg(data_policy)
  if ("data" %in% names(args_list) && !is.null(args_list$data) && policy != "none") {
    d <- args_list$data
    args_list$data <- switch(
      policy,
      names      = list(colnames = names(d), nrow = nrow(d)),
      signature  = utils::object.size(d),
      hash       = digest::digest(d, algo = algo, serialize = TRUE)
    )
  }

  digest::digest(args_list, algo = algo, serialize = TRUE)
}



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

hash_brms_call2 <- function(args_list, algo = "xxhash64") {
  requireNamespace("digest")
  digest::digest(keep, algo = algo)
}

hash_data <- function(data){
  if( is.data.frame(data ) && nrow(data) * ncol(data) > 1e7 ) {
    data_hash <- digest::digest(dim(data))
  }else{
    data_hash <- digest::digest(data )
  }
  data_hash
}

stancode_local <- function(...){

  dots <- nlist(... )

  other_dots <- dots
  other_dots$formula <- NULL
  other_dots$data <- NULL

  stancode(dots$formula, dots$data, rlang::splice(other_dots) )
}

stancode_local2 <- function(dots){

  other_dots <- dots
  other_dots$formula <- NULL
  other_dots$data <- NULL
  other_dots$family <- NULL
  other_dots$prior  <- NULL
  data <- dots$data
  family <- dots$family
  formula <- dots$formula
  data2 <-  dots$data2
  prior <- dots$prior
  sparse <-  dots$sparse
  autocor <-  dots$autocor
  cov_ranef <-  dots$cov_ranef
  drop_unused_levels <-  dots$drop_unused_levels
  sample_prior <-  dots$sample_prior

  stanvars <- dots$stanvars
  stan_funs <- dots$stan_funs

  .save_pars <- dots$save_pars
  save_ranef = dots$save_ranef
  save_mevars = dots$save_mevars
  save_all_pars = dots$save_all_pars
  normalize <- dots$normalize
  backend  <- dots$backend
  threads <- dots$threads
  save_model <- dots$save_model
  opencl <- dots$opencl

  formula <- validate_formula(
    formula, data = data, family = family,
    autocor = autocor, sparse = sparse,
    cov_ranef = cov_ranef
  )



  # family <- get_element(formula, "family")
  bterms <- brmsterms(formula)
  data2 <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor( formula ),
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
  .save_pars <- validate_save_pars(
    .save_pars, save_ranef = save_ranef,
    save_mevars = save_mevars,
    save_all_pars = save_all_pars
  )

  # generate Stan code
  model <- .stancode(
    bframe, prior = prior, stanvars = stanvars,
    save_model = save_model, backend = backend, threads = threads,
    opencl = opencl, normalize = normalize
  )



  model
  # stancode.default (dots$formula,  data = data  , family = family ,  rlang::splice(other_dots)  )


}


hash_brms_call_helper <-function(args_list){
  algo = "xxhash64"
  requireNamespace("digest")
  data <- args_list$data

  if( is.null(data) ){
    correct_data <-   epilepsy[  -c(1,8,11,25,29,12) , ]

  }else{
    correct_data <- data
  }



  args_list$data <- correct_data

  # code <- stancode(... , data = correct_data )
  # code <- stancode( !!!args_list , data = correct_data )
  # code <- do.call(stancode , args_list    )

  assign( "dbg_args_list" ,  args_list , envir = .GlobalEnv )

  code <- stancode_local2(args_list)

  formula <- as.character(args_list$formula)

  code <- as.character(code )



  brms_version <- packageVersion("brms")
  #handle data
  data  <- hash_data(correct_data)

  args_list <- nlist( code ,  formula , data ,  brms_version )

  digest::digest(args_list, algo = algo)



}

hash_brms_call <- function(...) {


  args_list <- brm_to_call(...)

  hash_brms_call_helper(args_list )
  # return(args_list)


}

..f <- function(){

  # hash_brms_call( count ~  zAge + zBase * Trt + (1|patient),
  #                 data = epilepsy, family = poisson() )
  # a = hash_brms_call( count ~  zAge + zBase * Trt + (1|patient),
  #                 data = epilepsy, family = poisson() )
  #
  # b = hash_brms_call( count ~  zAge + zBase * Trt + (1|patient),
  #                     data = epilepsy, family = poisson() )

  s1 <- function(formula  , ...){
    dots <- rlang::dots_list( ... )
    dots$formula <- formula

    rlang::splice(dots)
  }
  c1 <-   s1(count ~  zAge + zBase * Trt + (1|patient),
             data = epilepsy, family = poisson())

  brm_to_call(c1)

  a <-  brm_to_call( count ~  zAge + zBase * Trt + (1|patient),
                     data = epilepsy, family = poisson() )

  hash_brms_call_helper(a)

  # c1  <- rlang::dots_list(count ~  zAge + zBase * Trt + (1|patient),
  #         data = epilepsy, family = poisson())
  # brm_to_call( rlang::splice( c1  ))
  # hash_brms_call(rlang::splice( c1  ) )

  a <- hash_brms_call( count ~  zAge + zBase * Trt + (1|patient),
                       data = epilepsy, family = poisson() )

  b<- hash_brms_call( count ~  zBase * Trt + (1|patient),
                      data = epilepsy, family = poisson() )

  expect_false( a == b )

  a <- hash_brms_call( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  b <- hash_brms_call( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  expect_true( a == b )

  a <- hash_brms_call( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(2) , ], family = poisson() )

  b <- hash_brms_call( count ~  zBase * Trt + (1|patient),
                       data = epilepsy[-c(1) , ], family = poisson() )

  expect_false( a == b )



  m1 <- brm( count ~  zBase * Trt + (1|patient),
             data = epilepsy[-c(2) , ], family = poisson()  , file = "m")


}



# hash_dots
#
# Convenience wrapper for hash_brms_call that accepts ... instead of a named list.
# Useful when manually specifying formula, data, family, etc.
#
# Example:
#   hash_dots(formula = y ~ x, data = df, family = gaussian())
hash_dots  <- function(...){
  hash_brms_call( ...  )

  # dots <- nlist(...)
  # hash_brms_call( dots  )
}


hash_dots_general  <- function(...){


  dots <- nlist(...)
  hash_brms_call( dots  )
}





remove_env_attrs <- function(obj) {
  # Remove environment attribute if exists
  if (!is.null(attr(obj, ".Environment"))) {
    attr(obj, ".Environment") <- NULL
  }
  # If the object is a list or pairlist, recursively apply to elements
  if (is.list(obj) || is.pairlist(obj)) {
    obj <- lapply(obj, remove_env_attrs)
  }
  # If the object is a formula, also check its environment
  if (inherits(obj, "formula")) {
    environment(obj) <- NULL
  }
  obj
}

