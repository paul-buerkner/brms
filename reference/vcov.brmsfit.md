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
#> Chain 1: Gradient evaluation took 7.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.75 seconds.
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
#> Chain 1:  Elapsed Time: 0.639 seconds (Warm-up)
#> Chain 1:                0.695 seconds (Sampling)
#> Chain 1:                1.334 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.701 seconds (Warm-up)
#> Chain 2:                0.537 seconds (Sampling)
#> Chain 2:                1.238 seconds (Total)
#> Chain 2: 
#> Warning: There were 2 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
vcov(fit)
#>              Intercept         zAge        zBase        Trt1  zBase:Trt1
#> Intercept   0.88644211 -0.017835660  0.013806951 -0.59013560 -0.03719786
#> zAge       -0.01783566  0.261236750 -0.008156782  0.07865400  0.09320789
#> zBase       0.01380695 -0.008156782  0.592407870 -0.05244389 -0.61251781
#> Trt1       -0.59013560  0.078654001 -0.052443885  1.84160706  0.07998178
#> zBase:Trt1 -0.03719786  0.093207885 -0.612517806  0.07998178  1.12548974
# }
```
