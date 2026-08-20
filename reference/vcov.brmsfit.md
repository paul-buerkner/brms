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
#> Chain 1: Gradient evaluation took 3.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 1:  Elapsed Time: 0.672 seconds (Warm-up)
#> Chain 1:                0.494 seconds (Sampling)
#> Chain 1:                1.166 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.3e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.33 seconds.
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
#> Chain 2:  Elapsed Time: 0.687 seconds (Warm-up)
#> Chain 2:                0.387 seconds (Sampling)
#> Chain 2:                1.074 seconds (Total)
#> Chain 2: 
#> Warning: There were 14 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
vcov(fit)
#>              Intercept        zAge       zBase        Trt1  zBase:Trt1
#> Intercept   1.15032828 -0.04942021 -0.01158160 -0.80612017  0.01497913
#> zAge       -0.04942021  0.28527683 -0.01574603  0.12128221  0.11267743
#> zBase      -0.01158160 -0.01574603  0.54910225  0.05211377 -0.53902232
#> Trt1       -0.80612017  0.12128221  0.05211377  2.23626581  0.01782143
#> zBase:Trt1  0.01497913  0.11267743 -0.53902232  0.01782143  0.99071009
# }
```
