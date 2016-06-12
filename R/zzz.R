# Uncomment the code below to enable unit tests of S3 methods

brmsfit_example <- brm(count ~ Trt*Age + mono(Exp) + s(Age) +
                         offset(Age) + (1+Trt|visit),
                       data = data.frame(count = rpois(236, lambda = 20),
                                         visit = rep(1:4, each = 59),
                                         patient = rep(1:59, 4),
                                         Age = rnorm(236), Trt = rnorm(236),
                                         Exp = sample(1:5, 236, TRUE)),
                       family = gaussian(), sample_prior = TRUE,
                       autocor = cor_arma(~visit|patient, 1, 1),
                       prior = c(set_prior("normal(0,5)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       warmup = 10, iter = 40, chains = 2, testmode = TRUE)
brmsfit_example$fit@stanmodel <- new("stanmodel")


# Uncomment the code below to enable unit tests for new stan functions

# new_stan_functions <- function() {
#   # copy all new stan functions into a single .stan file and compile it
#   isystem <- system.file("chunks", package = "brms")
#   chunk_filenames <- list.files(isystem, pattern = "^fun_")
#   families <- list(cumulative("probit"), sratio("logit"),
#                    cratio("cloglog"), acat("cauchit"))
#   cse <- c(rep(FALSE, 2), rep(TRUE, 2))
#   ordinal_funs <- ulapply(seq_along(families), function(i)
#     stan_ordinal(families[[i]], cse = cse[i])$fun)
#   temp_file <- tempfile()
#   cat(paste0("functions { \n",
#              collapse("  #include '", chunk_filenames, "' \n"),
#              collapse(ordinal_funs), "} \nmodel {} \n"),
#       file = temp_file)
#   model <- rstan::stanc_builder(file = temp_file, isystem = isystem,
#                                 obfuscate_model_name = TRUE)
#   rstan::stan_model(stanc_ret = model)
# }
# new_stan_functions <- new_stan_functions()
