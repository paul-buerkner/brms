# Extract Variance and Correlation Components

This function calculates the estimated standard deviations, correlations
and covariances of the group-level terms in a multilevel model of class
`brmsfit`. For linear models, the residual standard deviations,
correlations and covariances are also returned.

## Usage

``` r
# S3 method for class 'brmsfit'
VarCorr(
  x,
  sigma = 1,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- x:

  An object of class `brmsfit`.

- sigma:

  Ignored (included for compatibility with
  [`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)).

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

  Currently ignored.

## Value

A list of lists (one per grouping factor), each with three elements: a
matrix containing the standard deviations, an array containing the
correlation matrix, and an array containing the covariance matrix with
variances on the diagonal.

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
#> Chain 1: Gradient evaluation took 6.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.61 seconds.
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
#> Chain 1:  Elapsed Time: 0.705 seconds (Warm-up)
#> Chain 1:                0.432 seconds (Sampling)
#> Chain 1:                1.137 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.735 seconds (Warm-up)
#> Chain 2:                0.423 seconds (Sampling)
#> Chain 2:                1.158 seconds (Total)
#> Chain 2: 
#> Warning: There were 44 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
VarCorr(fit)
#> $visit
#> $visit$sd
#>           Estimate Est.Error       Q2.5    Q97.5
#> Intercept 1.140915  1.094191 0.03764509 4.261134
#> Trt1      1.345855  1.266183 0.04438976 4.851442
#> 
#> $visit$cor
#> , , Intercept
#> 
#>             Estimate Est.Error       Q2.5     Q97.5
#> Intercept  1.0000000 0.0000000  1.0000000 1.0000000
#> Trt1      -0.1217731 0.6009979 -0.9743096 0.9333658
#> 
#> , , Trt1
#> 
#>             Estimate Est.Error       Q2.5     Q97.5
#> Intercept -0.1217731 0.6009979 -0.9743096 0.9333658
#> Trt1       1.0000000 0.0000000  1.0000000 1.0000000
#> 
#> 
#> $visit$cov
#> , , Intercept
#> 
#>             Estimate Est.Error          Q2.5     Q97.5
#> Intercept  2.4983440  4.938207   0.001417192 18.157265
#> Trt1      -0.7692375  3.509571 -11.125638255  2.243856
#> 
#> , , Trt1
#> 
#>             Estimate Est.Error          Q2.5     Q97.5
#> Intercept -0.7692375  3.509571 -11.125638255  2.243856
#> Trt1       3.4137435  7.023089   0.001970522 23.536494
#> 
#> 
#> 
#> $residual__
#> $residual__$sd
#>  Estimate Est.Error     Q2.5    Q97.5
#>  7.637655 0.3460858 7.007586 8.346337
#> 
#> 
# }
```
