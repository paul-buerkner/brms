#'
#' #' Collect brm() arguments into a clean list
#' #' @noRd
#' .brm_collect_args <- function(...) {
#'   call_env  <- parent.frame()
#'   arg_names <- names(formals(brm))
#'
#'   # Remove ... from arg_names because it's not bound and can't be accessed
#'   arg_names <- arg_names[arg_names != "..."]
#'
#'   arg_list <- setNames(
#'     lapply(arg_names, function(a) get(a, envir = call_env)),
#'     arg_names
#'   )
#'
#'   # Explicitly capture dots
#'   dot_args <- list(...)
#'   arg_list$dot_args <- dot_args
#'
#'   # Optional: validate here
#'   arg_list$algorithm <- match.arg(arg_list$algorithm, algorithm_choices())
#'   arg_list$backend   <- match.arg(arg_list$backend, backend_choices())
#'
#'   arg_list
#' }
#'
#'
#'
#' #' Bayesian Regression Models using 'Stan'
#' #' @export
#' brm <- function(formula, data= NULL, family = gaussian(), prior = NULL,
#'                 autocor = NULL, data2 = NULL, cov_ranef = NULL,
#'                 sample_prior = "no", sparse = NULL, knots = NULL,
#'                 drop_unused_levels = TRUE, stanvars = NULL, stan_funs = NULL,
#'                 fit = NA, save_pars = getOption("brms.save_pars", NULL),
#'                 save_ranef = NULL, save_mevars = NULL, save_all_pars = NULL,
#'                 init = NULL, inits = NULL, chains = 4,
#'                 iter = getOption("brms.iter", 2000),
#'                 warmup = floor(iter / 2), thin = 1,
#'                 cores = getOption("mc.cores", 1),
#'                 threads = getOption("brms.threads", NULL),
#'                 opencl = getOption("brms.opencl", NULL),
#'                 normalize = getOption("brms.normalize", TRUE),
#'                 control = NULL,
#'                 algorithm = getOption("brms.algorithm", "sampling"),
#'                 backend = getOption("brms.backend", "rstan"),
#'                 future = getOption("future", FALSE), silent = 1,
#'                 seed = NA, save_model = NULL, stan_model_args = list(),
#'                 file = NULL, file_compress = TRUE,
#'                 file_refit = getOption("brms.file_refit", "never"),
#'                 empty = FALSE, rename = TRUE, ...) {
#'
#'   args <- .brm_collect_args()
#'   # .brm_internal(args)
#' }
#'
#'
#' #' Check for an existing cached brmsfit and return it if valid
#' #' @noRd
#' .brm_check <- function(brm_call_list) {
#'   file       <- brm_call_list$file
#'   file_refit <- match.arg(brm_call_list$file_refit, file_refit_options())
#'
#'   # Load brmsfit only if refit is explicitly set to 'never'
#'   if (!is.null(file) && file_refit == "never") {
#'     x <- read_brmsfit(file)
#'     if (!is.null(x)) {
#'       return(x)
#'     }
#'   }
#'   NULL
#' }
#'
#'
#' #' Internal engine to evaluate and fit a brms model
#' #' @noRd
#' .brm_internal <- function(brm_call_list) {
#'
#'   # Check if fit object can be reused from file
#'   result <- .brm_check(brm_call_list)
#'   if (!is.null(result)) {
#'     return(result)
#'   }
#'
#'   # Extract parameters
#'   fit        <- brm_call_list$fit
#'   file       <- brm_call_list$file
#'   file_refit <- brm_call_list$file_refit
#'   algorithm  <- brm_call_list$algorithm
#'   backend    <- brm_call_list$backend
#'   silent     <- brm_call_list$silent
#'   data       <- brm_call_list$data
#'   data2      <- brm_call_list$data2
#'   family     <- brm_call_list$family
#'   model      <- brm_call_list$model
#'
#'   if (is.brmsfit(fit)) {
#'     # --- Reuse existing brmsfit object ---
#'     x <- fit
#'     x$criteria <- list()
#'     sdata <- standata(x)
#'
#'     if (!is.null(file) && file_refit == "on_change") {
#'       x_file <- read_brmsfit(file)
#'       if (!is.null(x_file)) {
#'         needs_refit <- brmsfit_needs_refit(
#'           x_file, scode = stancode(x), sdata = sdata,
#'           data = x$data, algorithm = algorithm, silent = silent
#'         )
#'         if (!needs_refit) {
#'           return(x_file)
#'         }
#'       }
#'     }
#'
#'     model   <- compiled_model(x)
#'     exclude <- exclude_pars(x)
#'
#'   } else {
#'     # --- Build a new brmsfit object from scratch ---
#'
#'     formula <- validate_formula(
#'       brm_call_list$formula, data, brm_call_list$family,
#'       brm_call_list$autocor, brm_call_list$sparse, brm_call_list$cov_ranef
#'     )
#'
#'     family <- get_element(formula, "family")
#'     bterms <- brmsterms(formula)
#'
#'     data2 <- validate_data2(
#'       brm_call_list$data2, bterms,
#'       get_data2_autocor(formula),
#'       get_data2_cov_ranef(formula)
#'     )
#'
#'     data <- validate_data(
#'       data, bterms, data2, brm_call_list$knots,
#'       brm_call_list$drop_unused_levels,
#'       data_name = substitute_name(data)
#'     )
#'
#'     bframe <- brmsframe(bterms, data)
#'
#'     prior <- .validate_prior(
#'       brm_call_list$prior, bframe, brm_call_list$sample_prior
#'     )
#'
#'     stanvars <- validate_stanvars(
#'       brm_call_list$stanvars, brm_call_list$stan_funs
#'     )
#'
#'     save_pars <- validate_save_pars(
#'       brm_call_list$save_pars,
#'       brm_call_list$save_ranef,
#'       brm_call_list$save_mevars,
#'       brm_call_list$save_all_pars
#'     )
#'
#'     model <- .stancode(
#'       bframe, prior, stanvars,
#'       brm_call_list$save_model, backend,
#'       brm_call_list$threads, brm_call_list$opencl,
#'       brm_call_list$normalize
#'     )
#'
#'     x <- brmsfit(
#'       formula = formula, data = data, data2 = data2,
#'       prior = prior, stanvars = stanvars, model = model,
#'       algorithm = algorithm, backend = backend,
#'       threads = brm_call_list$threads, opencl = brm_call_list$opencl,
#'       save_pars = save_pars, ranef = bframe$frame$re,
#'       family = family, basis = frame_basis(bframe, data),
#'       stan_args = nlist(
#'         init = brm_call_list$init,
#'         silent = silent,
#'         control = brm_call_list$control,
#'         stan_model_args = brm_call_list$stan_model_args,
#'         ... = brm_call_list$dot_args
#'       )
#'     )
#'
#'     exclude <- exclude_pars(x, bframe)
#'
#'     sdata <- .standata(
#'       bframe, data, prior, data2, stanvars,
#'       brm_call_list$threads
#'     )
#'
#'     if (brm_call_list$empty) return(x)
#'
#'     if (!is.null(file) && file_refit == "on_change") {
#'       x_file <- read_brmsfit(file)
#'       if (!is.null(x_file)) {
#'         needs_refit <- brmsfit_needs_refit(
#'           x_file, scode = model, sdata = sdata, data = data,
#'           algorithm = algorithm, silent = silent
#'         )
#'         if (!needs_refit) return(x_file)
#'       }
#'     }
#'
#'     # compile Stan model
#'     model <- do_call(compile_model, modifyList(
#'       brm_call_list$stan_model_args,
#'       list(model = model, backend = backend,
#'            threads = brm_call_list$threads,
#'            opencl = brm_call_list$opencl,
#'            silent = silent)
#'     ))
#'   }
#'
#'   # --- Fit Stan model ---
#'   fit_args <- nlist(
#'     model, sdata, algorithm, backend,
#'     iter     = brm_call_list$iter,
#'     warmup   = brm_call_list$warmup,
#'     thin     = brm_call_list$thin,
#'     chains   = brm_call_list$chains,
#'     cores    = brm_call_list$cores,
#'     threads  = brm_call_list$threads,
#'     opencl   = brm_call_list$opencl,
#'     init     = brm_call_list$init,
#'     exclude,
#'     control = brm_call_list$control,
#'     future   = brm_call_list$future,
#'     seed     = brm_call_list$seed,
#'     silent, ... = brm_call_list$...
#'   )
#'
#'   x$fit <- do_call(fit_model, fit_args)
#'
#'   if (brm_call_list$rename) {
#'     x <- rename_pars(x)
#'   }
#'
#'   if (!is.null(file)) {
#'     x <- write_brmsfit(x, file, compress = brm_call_list$file_compress)
#'   }
#'
#'   x
#' }
#'
