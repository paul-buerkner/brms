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
#> Chain 1: Gradient evaluation took 4.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.42 seconds.
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
#> Chain 1:  Elapsed Time: 0.642 seconds (Warm-up)
#> Chain 1:                0.39 seconds (Sampling)
#> Chain 1:                1.032 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
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
#> Chain 2:  Elapsed Time: 0.621 seconds (Warm-up)
#> Chain 2:                0.511 seconds (Sampling)
#> Chain 2:                1.132 seconds (Total)
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
#> 1  0.27600874 0.9067979 -1.296074 2.589330
#> 2  0.00516405 0.8384335 -1.747060 1.903566
#> 3  0.06650317 0.8422333 -1.689903 2.023606
#> 4 -0.22466639 0.8579908 -2.230082 1.435661
#> 
#> , , Trt1
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1  0.1836235  1.171169 -2.146618 3.089103
#> 2  0.2782745  1.156733 -1.803233 2.972164
#> 3  0.1664216  1.171482 -2.171578 2.982610
#> 4 -0.2708869  1.086226 -2.834679 1.930818
#> 
#> 
# }
```
