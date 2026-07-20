# Extract Group-Level Estimates

Extract the group-level ('random') effects of each level from a
`brmsfit` object.

## Usage

``` r
# S3 method for class 'brmsfit'
ranef(
  object,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  pars = NULL,
  groups = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

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

- pars:

  Optional names of coefficients to extract. By default, all
  coefficients are extracted.

- groups:

  Optional names of grouping variables for which to extract effects.

- ...:

  Currently ignored.

## Value

A list of 3D arrays (one per grouping factor). If `summary` is `TRUE`,
the 1st dimension contains the factor levels, the 2nd dimension contains
the summary statistics (see
[`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.md)),
and the 3rd dimension contains the group-level effects. If `summary` is
`FALSE`, the 1st dimension contains the posterior draws, the 2nd
dimension contains the factor levels, and the 3rd dimension contains the
group-level effects.

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zAge + zBase * Trt + (1+Trt|visit),
           data = epilepsy, family = gaussian(), chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.37 seconds.
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
#> Chain 1:  Elapsed Time: 0.727 seconds (Warm-up)
#> Chain 1:                0.485 seconds (Sampling)
#> Chain 1:                1.212 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.31 seconds.
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
#> Chain 2:  Elapsed Time: 0.678 seconds (Warm-up)
#> Chain 2:                0.656 seconds (Sampling)
#> Chain 2:                1.334 seconds (Total)
#> Chain 2: 
#> Warning: There were 4 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
ranef(fit)
#> $visit
#> , , Intercept
#> 
#>      Estimate Est.Error      Q2.5    Q97.5
#> 1  0.30467081 0.8750871 -1.214209 2.499196
#> 2  0.02825514 0.8187142 -1.690913 1.847562
#> 3  0.08762618 0.8016888 -1.578298 1.942982
#> 4 -0.21487766 0.8408538 -2.183174 1.372078
#> 
#> , , Trt1
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1  0.1678512  1.091097 -1.966271 2.887353
#> 2  0.2875135  1.123628 -1.756602 3.154275
#> 3  0.1174871  1.043511 -1.930729 2.721896
#> 4 -0.2282879  1.103589 -2.918115 1.886437
#> 
#> 
# }
```
