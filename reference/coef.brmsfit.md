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
#> Chain 1:  Elapsed Time: 0.722 seconds (Warm-up)
#> Chain 1:                0.629 seconds (Sampling)
#> Chain 1:                1.351 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.3 seconds.
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
#> Chain 2:  Elapsed Time: 0.663 seconds (Warm-up)
#> Chain 2:                0.481 seconds (Sampling)
#> Chain 2:                1.144 seconds (Total)
#> Chain 2: 
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
## extract population and group-level coefficients separately
fixef(fit)
#>              Estimate Est.Error       Q2.5     Q97.5
#> Intercept   8.4080560 0.9655297  6.3413222 10.130798
#> zAge        1.6382347 0.5237024  0.6278908  2.639194
#> zBase       7.0265912 0.7639634  5.4769891  8.482213
#> Trt1       -0.7099839 1.3194749 -3.3396086  1.723281
#> zBase:Trt1  4.7828378 1.0198244  2.8581711  6.749996
ranef(fit)
#> $visit
#> , , Intercept
#> 
#>      Estimate Est.Error      Q2.5    Q97.5
#> 1  0.34135610 0.9382166 -1.226836 2.772130
#> 2  0.09968401 0.8782179 -1.696261 2.183759
#> 3  0.14431020 0.9119715 -1.480372 2.288115
#> 4 -0.20093860 0.9044591 -2.163351 1.678969
#> 
#> , , Trt1
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1  0.2068644  1.129104 -1.896670 2.884655
#> 2  0.2626932  1.181336 -2.087958 3.195214
#> 3  0.1311467  1.152926 -2.244617 2.814154
#> 4 -0.2176979  1.209247 -3.083865 2.291449
#> 
#> 
## extract combined coefficients
coef(fit)
#> $visit
#> , , Intercept
#> 
#>   Estimate Est.Error     Q2.5     Q97.5
#> 1 8.749412 0.9369108 7.026388 10.724609
#> 2 8.507740 0.9026341 6.773999 10.419782
#> 3 8.552366 0.9153059 6.673310 10.425008
#> 4 8.207117 0.9397999 6.178842  9.968766
#> 
#> , , Trt1
#> 
#>     Estimate Est.Error      Q2.5    Q97.5
#> 1 -0.5031195  1.241760 -2.803109 2.036540
#> 2 -0.4472907  1.290644 -2.880706 2.143415
#> 3 -0.5788371  1.244170 -2.903376 1.944472
#> 4 -0.9276818  1.319342 -3.659297 1.526786
#> 
#> , , zAge
#> 
#>   Estimate Est.Error      Q2.5    Q97.5
#> 1 1.638235 0.5237024 0.6278908 2.639194
#> 2 1.638235 0.5237024 0.6278908 2.639194
#> 3 1.638235 0.5237024 0.6278908 2.639194
#> 4 1.638235 0.5237024 0.6278908 2.639194
#> 
#> , , zBase
#> 
#>   Estimate Est.Error     Q2.5    Q97.5
#> 1 7.026591 0.7639634 5.476989 8.482213
#> 2 7.026591 0.7639634 5.476989 8.482213
#> 3 7.026591 0.7639634 5.476989 8.482213
#> 4 7.026591 0.7639634 5.476989 8.482213
#> 
#> , , zBase:Trt1
#> 
#>   Estimate Est.Error     Q2.5    Q97.5
#> 1 4.782838  1.019824 2.858171 6.749996
#> 2 4.782838  1.019824 2.858171 6.749996
#> 3 4.782838  1.019824 2.858171 6.749996
#> 4 4.782838  1.019824 2.858171 6.749996
#> 
#> 
# }
```
