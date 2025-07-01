# Mock brm function to check file_auto parameter
brm_mock_for_file_auto <- function(formula, data = NULL , family = gaussian(), prior = NULL,
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

  # define file argument automatically when file_auto is TRUE
  if( file_auto ){
    orig_seed <- seed
    # This list must include only/all the parameters that may change the result
    args_list <- nlist(formula, data, family, prior, autocor, data2, cov_ranef,
                       sample_prior, sparse, knots, drop_unused_levels, stanvars,
                       stan_funs, fit, save_pars, save_ranef, save_mevars,
                       save_all_pars, init, inits, chains, iter, warmup, thin,
                       cores, threads, opencl, normalize, control, algorithm,
                       backend, future, orig_seed= orig_seed, stan_model_args, empty)
    # args_list <- match.call()
    auto_res <- create_filename_auto(file, file_refit, file_auto, args_list)
    file <- auto_res$file
    file_refit <- auto_res$file_refit
  }
  nlist(file, file_refit)
}
# file_auto 
test_that("file_auto option works", {
  # skip("Temporarily disabled for debugging reasons")
  # test_cache_dir <- tempdir()
  # options(brms.cache_folder = test_cache_dir)
    epilepsy2 <- epilepsy[-c(1), ]
    # same
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                           data = epilepsy, family = poisson() ,    file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = gaussian() ,    file_auto = TRUE)
    expect_equal(f1$file, f2$file)
    # different data
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = poisson(), file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy2, family = poisson(), file_auto = TRUE)
    expect_false(f1$file == f2$file)
    # different family
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = poisson(), file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = gaussian(), file_auto = TRUE)
    expect_equal(f1$file, f2$file)
})
