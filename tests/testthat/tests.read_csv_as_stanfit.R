context("Tests for read_csv_as_stanfit")

library(cmdstanr)

# fit model
model_code <- "
parameters {
  vector[2] x;
}
model {
  x ~ std_normal();
}
"

stan_file <- cmdstanr::write_stan_file(code = model_code)
mod <- cmdstan_model(stan_file)

fit <- mod$sample(parallel_chains = 2,
                  iter_warmup = 200,
                  iter_sampling=200)

fit_warmup <- mod$sample(parallel_chains = 2,
                         iter_warmup = 200,
                         iter_sampling=200,
                         save_warmup = T)

fit_dense_warmup <- mod$sample(parallel_chains = 4,
                               iter_warmup = 200,
                               iter_sampling=200,
                               metric = "dense_e",
                               save_warmup = T)

fit_nosampling <- mod$sample(parallel_chains = 4,
                             iter_warmup = 200,
                             iter_sampling = 0,
                             save_warmup = T)

fit_thinned <-  mod$sample(parallel_chains = 4,
                           iter_warmup = 200,
                           iter_sampling = 200,
                           thin = 5)

fit_variational <- mod$variational()

test_set <- list(
  single_chain = fit$output_files()[[1]],
  multi_chain = fit$output_files(),
  with_warmup = fit_warmup$output_files(),
  dense_warmup = fit_dense_warmup$output_files(),
  no_samples = fit_nosampling$output_files(),
  thinned = fit_thinned$output_files(),
  VI = fit_variational$output_files()
)

compare_functions <- function(filename, check_pars = TRUE) {
  rstan_read <- suppressWarnings(rstan::read_stan_csv(filename))
  csv_as_stanfit <- brms:::read_csv_as_stanfit(filename)


  # should only have permutation different so set to NULL
  rstan_read@sim$permutation <- NULL
  csv_as_stanfit@sim$permutation <- NULL

  if (check_pars) {
    # currently fails for VI because of different preprocessing
    expect_identical(rstan_read@model_pars, csv_as_stanfit@model_pars)
    expect_equal(rstan_read@par_dims, csv_as_stanfit@par_dims)
    expect_equal(rstan_read@sim, csv_as_stanfit@sim)
  }

  expect_identical(rstan_read@model_name, csv_as_stanfit@model_name)
  expect_identical(rstan_read@mode, csv_as_stanfit@mode)
  expect_equal(rstan_read@inits, csv_as_stanfit@inits)
  # should have 4 missing bits of info: metric_file, file, diagnostic_file, stancflags
  # expect_equal(length(all.equal(rstan_read@stan_args[[1]], csv_as_stanfit@stan_args[[1]])) == 4,
  expect_equal(rstan_read@stanmodel, csv_as_stanfit@stanmodel)
  expect_equal(rstan_read@date, csv_as_stanfit@date)
  return(invisible(NULL))
}

# tests
test_that("read methods identical: single chain with samples", {
  compare_functions(test_set$single_chain)
})

test_that("read methods identical: multiple chains with samples", {
  compare_functions(test_set$multi_chain)
})

test_that("read methods identical: warmup", {
  compare_functions(test_set$with_warmup)
})

test_that("read methods identical: dense warmup", {
  compare_functions(test_set$dense_warmup)
})

test_that("read methods identical: no samples", {
  compare_functions(test_set$no_samples)
})

test_that("read methods identical: thinned samples", {
  compare_functions(test_set$thinned)
})

test_that("read methods identical: variational inference", {
  # comparison of parameters and their draws may fail because
  # of CSV preprocessing done differently by rstan and cmdstanr
  compare_functions(test_set$VI, check_pars = FALSE)
})

