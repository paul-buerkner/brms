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
#> Chain 1: Gradient evaluation took 5.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.51 seconds.
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
#> Chain 1:  Elapsed Time: 2.063 seconds (Warm-up)
#> Chain 1:                2.036 seconds (Sampling)
#> Chain 1:                4.099 seconds (Total)
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
#> Chain 2:  Elapsed Time: 2.165 seconds (Warm-up)
#> Chain 2:                1.506 seconds (Sampling)
#> Chain 2:                3.671 seconds (Total)
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
#> Chain 3:  Elapsed Time: 2.18 seconds (Warm-up)
#> Chain 3:                1.568 seconds (Sampling)
#> Chain 3:                3.748 seconds (Total)
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
#> Chain 4:  Elapsed Time: 2.109 seconds (Warm-up)
#> Chain 4:                1.585 seconds (Sampling)
#> Chain 4:                3.694 seconds (Total)
#> Chain 4: 

# perform k-fold cross validation
(kf <- kfold(fit, save_fits = TRUE, chains = 1))
#> Fitting model 1 out of 10
#> Start sampling
#> Fitting model 2 out of 10
#> Start sampling
#> Fitting model 3 out of 10
#> Start sampling
#> Fitting model 4 out of 10
#> Start sampling
#> Fitting model 5 out of 10
#> Start sampling
#> Fitting model 6 out of 10
#> Start sampling
#> Fitting model 7 out of 10
#> Start sampling
#> Fitting model 8 out of 10
#> Start sampling
#> Fitting model 9 out of 10
#> Start sampling
#> Fitting model 10 out of 10
#> Start sampling
#> Warning: Found 4 observations with a pareto_k > 0.666666666666667 in model 'fit'.
#> 
#> Based on 10-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -653.2 33.8
#> p_kfold        75.5 13.5
#> kfoldic      1306.4 67.6
#> ------
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)     232   98.3%   290     
#>    (0.67, 1]   (bad)        3    1.3%   <NA>    
#>     (1, Inf)   (very bad)   1    0.4%   <NA>    
#> See help('pareto-k-diagnostic') for details.

# define a loss function
rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}

# predict responses and evaluate the loss
kfp <- kfold_predict(kf)
rmse(y = kfp$y, yrep = kfp$yrep)
#> [1] 6.240098
# }
```
