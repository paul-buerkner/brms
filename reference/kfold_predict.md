# Predictions from K-Fold Cross-Validation

Compute and evaluate predictions after performing K-fold
cross-validation via
[`kfold`](https://paulbuerkner.com/brms/reference/kfold.brmsfit.md).

## Usage

``` r
kfold_predict(x, method = "posterior_predict", resp = NULL, ...)
```

## Arguments

- x:

  Object of class `'kfold'` computed by
  [`kfold`](https://paulbuerkner.com/brms/reference/kfold.brmsfit.md).
  For `kfold_predict` to work, the fitted model objects need to have
  been stored via argument `save_fits` of
  [`kfold`](https://paulbuerkner.com/brms/reference/kfold.brmsfit.md).

- method:

  Method used to obtain predictions. Can be set to `"posterior_predict"`
  (the default), `"posterior_epred"`, or `"posterior_linpred"`. For more
  details, see the respective function documentations.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

A `list` with two slots named `'y'` and `'yrep'`. `y` is a named vector
of observed responses (names = row indices). `yrep` is an `array` of
predictions with the same structure as the chosen `method` would return
for all \\N\\ out-of-fold observations at once, with posterior draws on
the first dimension and observations on the second. See the
documentation of the underlying prediction function for details on
additional dimensions.

## See also

[`kfold`](https://paulbuerkner.com/brms/reference/kfold.brmsfit.md)

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zBase * Trt + (1|patient),
           data = epilepsy, family = poisson())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
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
#> Chain 1:  Elapsed Time: 2.087 seconds (Warm-up)
#> Chain 1:                2.057 seconds (Sampling)
#> Chain 1:                4.144 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.7e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.27 seconds.
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
#> Chain 2:  Elapsed Time: 2.192 seconds (Warm-up)
#> Chain 2:                1.519 seconds (Sampling)
#> Chain 2:                3.711 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.6e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
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
#> Chain 3:  Elapsed Time: 2.204 seconds (Warm-up)
#> Chain 3:                1.579 seconds (Sampling)
#> Chain 3:                3.783 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
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
#> Chain 4:  Elapsed Time: 2.127 seconds (Warm-up)
#> Chain 4:                1.601 seconds (Sampling)
#> Chain 4:                3.728 seconds (Total)
#> Chain 4: 

# perform k-fold cross validation
(kf <- kfold(fit, save_fits = TRUE, chains = 1))
#> Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 17 globals exported for future expression (‘FUN()’) is 617.46 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option The three largest globals are ‘FUN’ (219.54 MiB of class ‘function’), ‘up_args’ (201.23 MiB of class ‘list’) and ‘newdata’ (196.68 MiB of class ‘list’)

# define a loss function
rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}

# predict responses and evaluate the loss
kfp <- kfold_predict(kf)
#> Error: object 'kf' not found
rmse(y = kfp$y, yrep = kfp$yrep)
#> Error: object 'kfp' not found
# }
```
