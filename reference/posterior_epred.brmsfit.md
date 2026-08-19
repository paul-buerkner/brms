# Draws from the Expected Value of the Posterior Predictive Distribution

Compute posterior draws of the expected value of the posterior
predictive distribution. Can be performed for the data used to fit the
model (posterior predictive checks) or for new data. By definition,
these predictions have smaller variance than the posterior predictions
performed by the
[`posterior_predict.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md)
method. This is because only the uncertainty in the expected value of
the posterior predictive distribution is incorporated in the draws
computed by `posterior_epred` while the residual error is ignored there.
However, the estimated means of both methods averaged across draws
should be very similar.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_epred(
  object,
  newdata = NULL,
  re_formula = NULL,
  re.form = NULL,
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

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

  Optional name of a predicted distributional parameter. If specified,
  expected predictions of this parameters are returned.

- nlpar:

  Optional name of a predicted non-linear parameter. If specified,
  expected predictions of this parameters are returned.

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

## Value

An `array` of draws. For categorical and ordinal models, the output is
an S x N x C array. Otherwise, the output is an S x N matrix, where S is
the number of posterior draws, N is the number of observations, and C is
the number of categories. In multivariate models, an additional
dimension is added to the output which indexes along the different
response variables.

## Details

`NA` values within factors in `newdata`, are interpreted as if all dummy
variables of this factor are zero. This allows, for instance, to make
predictions of the grand mean when using sum coding.

In multilevel models, it is possible to allow new levels of grouping
factors to be used in the predictions. This can be controlled via
argument `allow_new_levels`. New levels can be sampled in multiple ways,
which can be controlled via argument `sample_new_levels`. Both of these
arguments are documented in
[`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
along with several other useful arguments to control specific aspects of
the predictions.

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
#> Chain 1: Gradient evaluation took 6.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.68 seconds.
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
#> Chain 1:  Elapsed Time: 1.548 seconds (Warm-up)
#> Chain 1:                0.776 seconds (Sampling)
#> Chain 1:                2.324 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.44 seconds.
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
#> Chain 2:  Elapsed Time: 1.561 seconds (Warm-up)
#> Chain 2:                0.775 seconds (Sampling)
#> Chain 2:                2.336 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 7.2e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.72 seconds.
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
#> Chain 3:  Elapsed Time: 1.64 seconds (Warm-up)
#> Chain 3:                0.776 seconds (Sampling)
#> Chain 3:                2.416 seconds (Total)
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
#> Chain 4:  Elapsed Time: 1.599 seconds (Warm-up)
#> Chain 4:                0.773 seconds (Sampling)
#> Chain 4:                2.372 seconds (Total)
#> Chain 4: 

## compute expected predictions
ppe <- posterior_epred(fit)
str(ppe)
#>  num [1:4000, 1:572] 1.566 0.992 1.419 1.386 1.368 ...
# }
```
