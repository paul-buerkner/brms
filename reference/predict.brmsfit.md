# Draws from the Posterior Predictive Distribution

This method is an alias of
[`posterior_predict.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md)
with additional arguments for obtaining summaries of the computed draws.

## Usage

``` r
# S3 method for class 'brmsfit'
predict(
  object,
  newdata = NULL,
  re_formula = NULL,
  transform = NULL,
  resp = NULL,
  negative_rt = FALSE,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ntrys = 5,
  cores = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
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

- transform:

  (Deprecated) A function or a character string naming a function to be
  applied on the predicted responses before summary statistics are
  computed.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- negative_rt:

  Only relevant for Wiener diffusion models. A flag indicating whether
  response times of responses on the lower boundary should be returned
  as negative values. This allows to distinguish responses on the upper
  and lower boundary. Defaults to `FALSE`.

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

- ntrys:

  Parameter used in rejection sampling for truncated discrete models
  only (defaults to `5`). See Details for more information.

- cores:

  Number of cores (defaults to `1`). On non-Windows systems, this
  argument can be set globally via the `mc.cores` option.

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`.

- probs:

  The percentiles to be computed by the `quantile` function. Only used
  if `summary` is `TRUE`.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

An `array` of predicted response values. If `summary = FALSE` the output
resembles those of
[`posterior_predict.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md).

If `summary = TRUE` the output depends on the family: For categorical
and ordinal families, the output is an N x C matrix, where N is the
number of observations, C is the number of categories, and the values
are predicted category probabilities. For all other families, the output
is a N x E matrix where E = `2 + length(probs)` is the number of summary
statistics: The `Estimate` column contains point estimates (either mean
or median depending on argument `robust`), while the `Est.Error` column
contains uncertainty estimates (either standard deviation or median
absolute deviation depending on argument `robust`). The remaining
columns starting with `Q` contain quantile estimates as specified via
argument `probs`.

## See also

[`posterior_predict.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md)

## Examples

``` r
# \dontrun{
## fit a model
fit <- brm(time | cens(censored) ~ age + sex + (1 + age || patient),
           data = kidney, family = "exponential", init = "0")
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.33 seconds.
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
#> Chain 1:  Elapsed Time: 1.387 seconds (Warm-up)
#> Chain 1:                0.661 seconds (Sampling)
#> Chain 1:                2.048 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
#> Chain 2:  Elapsed Time: 1.271 seconds (Warm-up)
#> Chain 2:                0.661 seconds (Sampling)
#> Chain 2:                1.932 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 3:  Elapsed Time: 1.383 seconds (Warm-up)
#> Chain 3:                0.662 seconds (Sampling)
#> Chain 3:                2.045 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.2e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
#> Chain 4:  Elapsed Time: 1.368 seconds (Warm-up)
#> Chain 4:                0.66 seconds (Sampling)
#> Chain 4:                2.028 seconds (Total)
#> Chain 4: 

## predicted responses
pp <- predict(fit)
head(pp)
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,]  35.34280  47.65143 0.6382854  170.6624
#> [2,] 160.42088 296.89426 2.4560518  802.2983
#> [3,]  41.34877  55.35052 0.7301142  186.6020
#> [4,] 296.71247 384.85454 5.5359669 1221.5351
#> [5,]  46.87249  68.61577 0.8644541  225.8016
#> [6,] 198.86530 255.89102 3.7154450  910.0475

## predicted responses excluding the group-level effect of age
pp <- predict(fit, re_formula = ~ (1 | patient))
head(pp)
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,]  39.90130  52.84751 0.6616917  182.3723
#> [2,] 170.28120 595.24640 3.3010333  722.7105
#> [3,]  43.76757  56.27128 0.9051928  192.3669
#> [4,] 259.75823 324.36563 4.6723401 1087.5806
#> [5,]  48.62304  79.00334 0.9731580  242.6778
#> [6,] 197.38136 253.26588 3.4668561  878.3151

## predicted responses of patient 1 for new data
newdata <- data.frame(
  sex = factor(c("male", "female")),
  age = c(20, 50),
  patient = c(1, 1)
)
predict(fit, newdata = newdata)
#>       Estimate Est.Error      Q2.5    Q97.5
#> [1,]  37.68225  51.47541 0.7504945 164.7520
#> [2,] 142.02674 248.87428 1.9691851 711.8125
# }
```
