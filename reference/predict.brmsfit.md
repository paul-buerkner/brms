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
  output = "random",
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

- output:

  Type of predictive quantity to return. Either `"random"` (the
  default), `"probability"`, `"pit"`, `"density"`, or `"quantile"`. See
  Details for more information.

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
#> Chain 1: Gradient evaluation took 3.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.32 seconds.
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
#> Chain 1:  Elapsed Time: 1.126 seconds (Warm-up)
#> Chain 1:                0.582 seconds (Sampling)
#> Chain 1:                1.708 seconds (Total)
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
#> Chain 2:  Elapsed Time: 1.232 seconds (Warm-up)
#> Chain 2:                0.584 seconds (Sampling)
#> Chain 2:                1.816 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.2e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
#> Chain 3:  Elapsed Time: 1.193 seconds (Warm-up)
#> Chain 3:                0.578 seconds (Sampling)
#> Chain 3:                1.771 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
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
#> Chain 4:  Elapsed Time: 1.221 seconds (Warm-up)
#> Chain 4:                0.585 seconds (Sampling)
#> Chain 4:                1.806 seconds (Total)
#> Chain 4: 

## predicted responses
pp <- predict(fit)
head(pp)
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,]  34.83945  45.08492 0.6235013  159.3236
#> [2,] 153.64687 226.89030 2.6384661  715.0668
#> [3,]  40.52708  51.94876 0.7590829  187.9565
#> [4,] 292.03112 367.65340 5.7617542 1224.3354
#> [5,]  45.34998  62.82212 0.8670641  205.5986
#> [6,] 197.71048 262.87404 3.5950602  881.8410

## predicted responses excluding the group-level effect of age
pp <- predict(fit, re_formula = ~ (1 | patient))
head(pp)
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,]  39.36460  49.63430 0.6788658  176.2415
#> [2,] 165.96116 215.54506 3.0799573  761.5725
#> [3,]  43.03543  53.56215 0.8635685  181.4419
#> [4,] 256.42133 316.53610 5.0372193 1105.5024
#> [5,]  46.98695  64.98686 0.9155307  221.9494
#> [6,] 199.63181 262.03207 3.5420339  910.0072

## predicted responses of patient 1 for new data
newdata <- data.frame(
  sex = factor(c("male", "female")),
  age = c(20, 50),
  patient = c(1, 1)
)
predict(fit, newdata = newdata)
#>       Estimate Est.Error      Q2.5    Q97.5
#> [1,]  38.43111  54.22323 0.6979168 173.2195
#> [2,] 138.92795 210.06104 1.5201422 676.9411
# }
```
