# Draws from the Posterior Predictive Distribution

Compute posterior draws of the posterior predictive distribution. Can be
performed for the data used to fit the model (posterior predictive
checks) or for new data. By definition, these draws have higher variance
than draws of the expected value of the posterior predictive
distribution computed by
[`posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md).
This is because the residual error is incorporated in
`posterior_predict`. However, the estimated means of both methods
averaged across draws should be very similar.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_predict(
  object,
  newdata = NULL,
  re_formula = NULL,
  re.form = NULL,
  transform = NULL,
  resp = NULL,
  negative_rt = FALSE,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ntrys = 5,
  output = "random",
  q = NULL,
  p = NULL,
  cores = NULL,
  lower.tail = TRUE,
  log.p = FALSE,
  log = FALSE,
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

- output:

  Type of predictive quantity to return. Either `"random"` (the
  default), `"probability"`, `"pit"`, `"density"`, or `"quantile"`. See
  Details for more information.

- q:

  Optional values at which to evaluate probabilities, PIT values, or
  densities. Only used if `output` is `"probability"`, `"pit"`, or
  `"density"`. Defaults to the observed responses. See Details.

- p:

  Probabilities at which to evaluate quantiles. Only used and required
  if `output = "quantile"`.

- cores:

  Number of cores (defaults to `1`). On non-Windows systems, this
  argument can be set globally via the `mc.cores` option.

- lower.tail:

  Logical; If `TRUE` (default), return \\P(X \le x)\\ or the
  corresponding lower-tail quantiles. Else, return \\P(X \> x)\\ or
  upper-tail quantiles. Only used if `output` is `"probability"`,
  `"pit"`, or `"quantile"`.

- log.p:

  Logical; If `TRUE`, probabilities are given and returned on the log
  scale. Only used if `output` is `"probability"`, `"pit"`, or
  `"quantile"`.

- log:

  Logical; If `TRUE`, densities are returned on the log scale. Only used
  if `output = "density"`.

- ...:

  Further arguments passed to
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.md)
  that control several aspects of data validation and prediction.

## Value

An `array` of posterior predictive draws or related quantities (see
argument `output`). In univariate models, the output is an S x N matrix,
where S is the number of posterior draws and N is the number of
observations. In multivariate models, an additional dimension is added
to the output which indexes along the different response variables.

## Details

Besides random draws, related quantities such as CDF or PIT values,
densities, or quantiles can be returned via argument `output`. This
argument controls which predictive quantity is returned for each
posterior draw and observation:

- `"random"`: Draws from the posterior predictive distribution. This is
  the default and reproduces the historical behavior of
  `posterior_predict`.

- `"probability"`: Values of the posterior predictive CDF evaluated at
  `q`.

- `"pit"`: Probability integral transform (PIT) values evaluated at `q`.
  For continuous distributions, this is identical to `"probability"`.
  For discrete distributions, a randomized PIT is used to avoid point
  masses at the CDF jumps (Czado et al., 2009). See Säilynoja et
  al. (2022) and Tesso and Vehtari (2026) for using PIT values in
  predictive model checking.

- `"density"`: Values of the posterior predictive density or probability
  mass function evaluated at `q`.

- `"quantile"`: Quantiles of the posterior predictive distribution at
  probabilities `p`.

If `q` is `NULL`, the observed response values are used. Argument `p`
must be supplied when `output = "quantile"`. Not all response families
support outputs other than `"random"` yet; requesting an unsupported
combination raises an error.

For truncated discrete models only: In the absence of any general
algorithm to sample from truncated discrete distributions, rejection
sampling is applied in this special case. This means that values are
sampled until a value lies within the defined truncation boundaries. In
practice, this procedure may be rather slow (especially in R). Thus, we
try to do approximate rejection sampling by sampling each value `ntrys`
times and then select a valid value. If all values are invalid, the
closest boundary is used, instead. If there are more than a few of these
pathological cases, a warning will occur suggesting to increase argument
`ntrys`.

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

## References

Czado, C., Gneiting, T., & Held, L. (2009). Predictive model assessment
for count data. *Biometrics*, 65(4), 1254-1261.
[doi:10.1111/j.1541-0420.2009.01191.x](https://doi.org/10.1111/j.1541-0420.2009.01191.x)

Säilynoja, T., Bürkner, P.-C., & Vehtari, A. (2022). Graphical test for
discrete uniformity and its applications in goodness-of-fit evaluation
and multiple sample comparison. *Statistics and Computing*, 32, 32.
[doi:10.1007/s11222-022-10090-6](https://doi.org/10.1007/s11222-022-10090-6)

Tesso, H., & Vehtari, A. (2026). LOO-PIT predictive model checking.
*arXiv preprint*.
[doi:10.48550/arXiv.2603.02928](https://doi.org/10.48550/arXiv.2603.02928)

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
#> Chain 1: Gradient evaluation took 4.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
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
#> Chain 1:  Elapsed Time: 1.283 seconds (Warm-up)
#> Chain 1:                0.613 seconds (Sampling)
#> Chain 1:                1.896 seconds (Total)
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
#> Chain 2:  Elapsed Time: 1.271 seconds (Warm-up)
#> Chain 2:                0.444 seconds (Sampling)
#> Chain 2:                1.715 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.8e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
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
#> Chain 3:  Elapsed Time: 1.257 seconds (Warm-up)
#> Chain 3:                0.607 seconds (Sampling)
#> Chain 3:                1.864 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.8e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
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
#> Chain 4:  Elapsed Time: 1.403 seconds (Warm-up)
#> Chain 4:                0.617 seconds (Sampling)
#> Chain 4:                2.02 seconds (Total)
#> Chain 4: 

## predicted responses
pp <- posterior_predict(fit)
str(pp)
#>  num [1:4000, 1:76] 161.24 4.56 20.58 2.71 20.54 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL

## predicted responses excluding the group-level effect of age
pp <- posterior_predict(fit, re_formula = ~ (1 | patient))
str(pp)
#>  num [1:4000, 1:76] 93.79 14.89 34.25 19.37 3.91 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL

## predicted responses of patient 1 for new data
newdata <- data.frame(
  sex = factor(c("male", "female")),
  age = c(20, 50),
  patient = c(1, 1)
)
pp <- posterior_predict(fit, newdata = newdata)
str(pp)
#>  num [1:4000, 1:2] 315.39 7.98 31.82 10.96 73.23 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL
# }
```
