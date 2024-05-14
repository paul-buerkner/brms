# Uncomment the code below to enable unit tests for new stan functions

# new_stan_functions <- function() {
#   # copy all new stan functions into a single .stan file and compile it
#   isystem <- system.file("chunks", package = "brms")
#   chunk_filenames <- list.files(isystem, pattern = "^fun_")
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
