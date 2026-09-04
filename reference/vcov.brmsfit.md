# Covariance and Correlation Matrix of Population-Level Effects

Get a point estimate of the covariance or correlation matrix of
population-level parameters

## Usage

``` r
# S3 method for class 'brmsfit'
vcov(object, correlation = FALSE, pars = NULL, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- correlation:

  Logical; if `FALSE` (the default), compute the covariance matrix, if
  `TRUE`, compute the correlation matrix.

- pars:

  Optional names of coefficients to extract. By default, all
  coefficients are extracted.

- ...:

  Currently ignored.

## Value

covariance or correlation matrix of population-level parameters

## Details

Estimates are obtained by calculating the maximum likelihood covariances
(correlations) of the posterior draws.

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
#> Chain 1:  Elapsed Time: 0.663 seconds (Warm-up)
#> Chain 1:                0.584 seconds (Sampling)
#> Chain 1:                1.247 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 6.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.64 seconds.
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
#> Chain 2:  Elapsed Time: 0.65 seconds (Warm-up)
#> Chain 2:                0.396 seconds (Sampling)
#> Chain 2:                1.046 seconds (Total)
#> Chain 2: 
#> Warning: There were 7 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
vcov(fit)
#>              Intercept         zAge       zBase         Trt1  zBase:Trt1
#> Intercept   0.93935805 -0.011653242 -0.02894459 -0.536356833  0.03132944
#> zAge       -0.01165324  0.259764792 -0.03222474  0.005324015  0.13069852
#> zBase      -0.02894459 -0.032224743  0.61246274 -0.039686536 -0.64777809
#> Trt1       -0.53635683  0.005324015 -0.03968654  1.914384150  0.04634288
#> zBase:Trt1  0.03132944  0.130698518 -0.64777809  0.046342881  1.18649881
# }
```
