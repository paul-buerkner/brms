# Uncomment the code below to enable automated unit tests of S3 methods

brmsfit_example <- brm(count ~ Trt_c + offset(log_Age_c) + (1+Trt_c|visit),
                       data = data.frame(count = rpois(236, lambda = 20),
                                         visit = rep(1:4, each = 59),
                                         patient = rep(1:59, 4),
                                         log_Age_c = rnorm(236),
                                         Trt_c = rnorm(236)), 
                       family = gaussian(), sample_prior = TRUE,
                       autocor = cor_arma(~visit|patient, 1, 1),
                       prior = c(set_prior("normal(0,5)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       warmup = 10, iter = 40, chains = 2, testmode = TRUE)
brmsfit_example$fit@stanmodel <- new("stanmodel")


# Uncomment the code below to enable unit tests for new stan functions

new_stan_functions <- function() {
  # copy all new stan functions into a single .stan file and compile it 
  chunk_filenames <- list.files(system.file("chunks", package = "brms"))
  ordinal_funs <- ulapply(list(cumulative(), sratio(), cratio(), acat()),
    function(fam) stan_ordinal(fam, partial = TRUE)$fun)
  temp_file <- tempfile()
  cat(paste0("functions { \n",
             collapse("  #include '", chunk_filenames, "' \n"),
             collapse(ordinal_funs), "} \nmodel {} \n"), 
      file = temp_file)
  isystem <- system.file("chunks", package = "brms")
  model <- rstan::stanc_builder(file = temp_file, isystem = isystem,
                                obfuscate_model_name = TRUE)
  rstan::stan_model(stanc_ret = model)
}
new_stan_functions <- new_stan_functions()
