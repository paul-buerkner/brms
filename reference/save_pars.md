# Control Saving of Parameter Draws

Control which (draws of) parameters should be saved in a brms model. The
output of this function is meant for usage in the `save_pars` argument
of [`brm`](https://paulbuerkner.com/brms/reference/brm.md).

## Usage

``` r
save_pars(group = TRUE, latent = FALSE, all = FALSE, manual = NULL)
```

## Arguments

- group:

  A flag to indicate if group-level coefficients for each level of the
  grouping factors should be saved (default is `TRUE`). Set to `FALSE`
  to save memory. Alternatively, `group` may also be a character vector
  naming the grouping factors for which to save draws of coefficients.

- latent:

  A flag to indicate if draws of latent variables obtained by using `me`
  and `mi` terms should be saved (default is `FALSE`). Saving these
  draws allows to better use methods such as `posterior_predict` with
  the latent variables but leads to very large R objects even for models
  of moderate size and complexity. Alternatively, `latent` may also be a
  character vector naming the latent variables for which to save draws.

- all:

  A flag to indicate if draws of all variables defined in Stan's
  `parameters` block should be saved (default is `FALSE`). Saving these
  draws is required in order to apply the certain methods such as
  `bridge_sampler` and `bayes_factor`.

- manual:

  A character vector naming Stan variable names which should be saved.
  These names should match the variable names inside the Stan code
  before renaming. This feature is meant for power users only and will
  rarely be useful outside of very special cases.

## Value

A list of class `"save_pars"`.

## Examples

``` r
# \dontrun{
# don't store group-level coefficients
fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
           data = epilepsy, family = poisson(),
           save_pars = save_pars(group = FALSE))
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.6 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.985 seconds (Warm-up)
#> Chain 1:                1.786 seconds (Sampling)
#> Chain 1:                3.771 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.6e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.879 seconds (Warm-up)
#> Chain 2:                1.534 seconds (Sampling)
#> Chain 2:                3.413 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.4e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.24 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 2.101 seconds (Warm-up)
#> Chain 3:                1.414 seconds (Sampling)
#> Chain 3:                3.515 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.7e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.27 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 2.059 seconds (Warm-up)
#> Chain 4:                1.415 seconds (Sampling)
#> Chain 4:                3.474 seconds (Total)
#> Chain 4: 
variables(fit)
#> [1] "b_Intercept"           "b_zAge"                "b_zBase"              
#> [4] "b_Trt1"                "b_zBase:Trt1"          "sd_patient__Intercept"
#> [7] "Intercept"             "lprior"                "lp__"                 
# }
```
