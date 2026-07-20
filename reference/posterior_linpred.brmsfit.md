# Posterior Draws of the Linear Predictor

Compute posterior draws of the linear predictor, that is draws before
applying any link functions or other transformations. Can be performed
for the data used to fit the model (posterior predictive checks) or for
new data.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  re_formula = NULL,
  re.form = NULL,
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  incl_thres = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- transform:

  Logical; if `FALSE` (the default), draws of the linear predictor are
  returned. If `TRUE`, draws of the transformed linear predictor, that
  is, after applying the inverse link function are returned.

- newdata:

  An optional data.frame for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. `NA` values within
  factors (excluding grouping variables) are interpreted as if all dummy
  variables of this factor are zero. This allows, for instance, to make
  predictions of the grand mean when using sum coding. `NA` values
  within grouping variables are treated as a new level.

- re_formula:

  formula containing group-level effects to be considered in the
  prediction. If `NULL` (default), include all group-level effects; if
  `NA` or `~0`, include no group-level effects.

- re.form:

  Alias of `re_formula`.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- dpar:

  Name of a predicted distributional parameter for which draws are to be
  returned. By default, draws of the main distributional parameter(s)
  `"mu"` are returned.

- nlpar:

  Optional name of a predicted non-linear parameter. If specified,
  expected predictions of this parameters are returned.

- incl_thres:

  Logical; only relevant for ordinal models when `transform` is `FALSE`,
  and ignored otherwise. Shall the thresholds and category-specific
  effects be included in the linear predictor? For backwards
  compatibility, the default is to not include them.

- ndraws:

  Positive integer indicating how many posterior draws should be used.
  If `NULL` (the default) all draws are used. Ignored if `draw_ids` is
  not `NULL`.

- draw_ids:

  An integer vector specifying the posterior draws to be used. If `NULL`
  (the default), all draws are used.

- sort:

  Logical. Only relevant for time series models. Indicating whether to
  return predicted values in the original order (`FALSE`; default) or in
  the order of the time series (`TRUE`).

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## See also

[`posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md)

## Examples

``` r
# \dontrun{
## fit a model
fit <- brm(rating ~ treat + period + carry + (1|subject),
           data = inhaler)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.52 seconds.
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
#> Chain 1:  Elapsed Time: 1.617 seconds (Warm-up)
#> Chain 1:                0.771 seconds (Sampling)
#> Chain 1:                2.388 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 2:  Elapsed Time: 1.555 seconds (Warm-up)
#> Chain 2:                0.767 seconds (Sampling)
#> Chain 2:                2.322 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 3:  Elapsed Time: 1.539 seconds (Warm-up)
#> Chain 3:                0.76 seconds (Sampling)
#> Chain 3:                2.299 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4.5e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 4:  Elapsed Time: 1.543 seconds (Warm-up)
#> Chain 4:                0.757 seconds (Sampling)
#> Chain 4:                2.3 seconds (Total)
#> Chain 4: 
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

## extract linear predictor values
pl <- posterior_linpred(fit)
str(pl)
#>  num [1:4000, 1:572] 1.459 0.968 1.187 1.197 1.108 ...
# }
```
