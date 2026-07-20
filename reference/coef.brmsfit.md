# Extract Model Coefficients

Extract model coefficients, which are the sum of population-level
effects and corresponding group-level effects

## Usage

``` r
# S3 method for class 'brmsfit'
coef(object, summary = TRUE, robust = FALSE, probs = c(0.025, 0.975), ...)
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

- ...:

  Further arguments passed to
  [`fixef.brmsfit`](https://paulbuerkner.com/brms/reference/fixef.brmsfit.md)
  and
  [`ranef.brmsfit`](https://paulbuerkner.com/brms/reference/ranef.brmsfit.md).

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
#> Chain 1: Gradient evaluation took 4.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.43 seconds.
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
#> Chain 1:  Elapsed Time: 0.747 seconds (Warm-up)
#> Chain 1:                0.398 seconds (Sampling)
#> Chain 1:                1.145 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.6e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.36 seconds.
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
#> Chain 2:  Elapsed Time: 0.756 seconds (Warm-up)
#> Chain 2:                0.57 seconds (Sampling)
#> Chain 2:                1.326 seconds (Total)
#> Chain 2: 
#> Warning: There were 8 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
## extract population and group-level coefficients separately
fixef(fit)
#>              Estimate Est.Error       Q2.5     Q97.5
#> Intercept   8.3805271 1.0295451  6.1895980 10.415388
#> zAge        1.6255192 0.5033540  0.6195475  2.600940
#> zBase       7.0736887 0.7652248  5.5799384  8.553424
#> Trt1       -0.5699312 1.3036311 -3.2956938  2.031424
#> zBase:Trt1  4.7561967 1.0215380  2.7296836  6.723061
ranef(fit)
#> $visit
#> , , Intercept
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1  0.3285308 0.9344019 -1.293782 2.693919
#> 2  0.0750183 0.9101700 -1.837578 2.313599
#> 3  0.1297632 0.8784518 -1.587116 2.263351
#> 4 -0.2364204 0.9471601 -2.552662 1.704648
#> 
#> , , Trt1
#> 
#>      Estimate Est.Error      Q2.5    Q97.5
#> 1  0.13890322  1.118502 -2.098429 2.660623
#> 2  0.21155455  1.128022 -1.794073 3.029122
#> 3  0.08166616  1.058476 -2.116232 2.612850
#> 4 -0.31902731  1.137901 -3.214865 1.670180
#> 
#> 
## extract combined coefficients
coef(fit)
#> $visit
#> , , Intercept
#> 
#>   Estimate Est.Error     Q2.5    Q97.5
#> 1 8.709058 0.9577679 6.942009 10.71741
#> 2 8.455545 0.9243263 6.648054 10.24409
#> 3 8.510290 0.9730510 6.580128 10.47897
#> 4 8.144107 0.9591460 6.149919 10.06692
#> 
#> , , Trt1
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1 -0.4310280  1.293600 -3.089817 2.101123
#> 2 -0.3583766  1.296555 -2.783107 2.205693
#> 3 -0.4882650  1.281349 -3.041544 1.993971
#> 4 -0.8889585  1.334815 -3.735598 1.574951
#> 
#> , , zAge
#> 
#>   Estimate Est.Error      Q2.5   Q97.5
#> 1 1.625519  0.503354 0.6195475 2.60094
#> 2 1.625519  0.503354 0.6195475 2.60094
#> 3 1.625519  0.503354 0.6195475 2.60094
#> 4 1.625519  0.503354 0.6195475 2.60094
#> 
#> , , zBase
#> 
#>   Estimate Est.Error     Q2.5    Q97.5
#> 1 7.073689 0.7652248 5.579938 8.553424
#> 2 7.073689 0.7652248 5.579938 8.553424
#> 3 7.073689 0.7652248 5.579938 8.553424
#> 4 7.073689 0.7652248 5.579938 8.553424
#> 
#> , , zBase:Trt1
#> 
#>   Estimate Est.Error     Q2.5    Q97.5
#> 1 4.756197  1.021538 2.729684 6.723061
#> 2 4.756197  1.021538 2.729684 6.723061
#> 3 4.756197  1.021538 2.729684 6.723061
#> 4 4.756197  1.021538 2.729684 6.723061
#> 
#> 
# }
```
