# Compute posterior uncertainty intervals

Compute posterior uncertainty intervals for `brmsfit` objects.

## Usage

``` r
# S3 method for class 'brmsfit'
posterior_interval(object, pars = NA, variable = NULL, prob = 0.95, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- pars:

  Deprecated alias of `variable`. For reasons of backwards
  compatibility, `pars` is interpreted as a vector of regular
  expressions by default unless `fixed = TRUE` is specified.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- prob:

  A value between 0 and 1 indicating the desired probability to be
  covered by the uncertainty intervals. The default is 0.95.

- ...:

  More arguments passed to
  [`as.matrix.brmsfit`](https://paulbuerkner.com/brms/reference/as.data.frame.brmsfit.md).

## Value

A `matrix` with lower and upper interval bounds as columns and as many
rows as selected variables.

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zAge + zBase * Trt,
           data = epilepsy, family = negbinomial())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
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
#> Chain 1:  Elapsed Time: 0.179 seconds (Warm-up)
#> Chain 1:                0.162 seconds (Sampling)
#> Chain 1:                0.341 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.5e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.25 seconds.
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
#> Chain 2:  Elapsed Time: 0.181 seconds (Warm-up)
#> Chain 2:                0.148 seconds (Sampling)
#> Chain 2:                0.329 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.25 seconds.
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
#> Chain 3:  Elapsed Time: 0.182 seconds (Warm-up)
#> Chain 3:                0.151 seconds (Sampling)
#> Chain 3:                0.333 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.3e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 4:  Elapsed Time: 0.171 seconds (Warm-up)
#> Chain 4:                0.164 seconds (Sampling)
#> Chain 4:                0.335 seconds (Total)
#> Chain 4: 
posterior_interval(fit)
#>                       2.5%         97.5%
#> b_Intercept   1.751540e+00  2.038595e+00
#> b_zAge        4.318639e-03  2.186405e-01
#> b_zBase       5.753462e-01  8.977808e-01
#> b_Trt1       -3.883649e-01  7.795978e-03
#> b_zBase:Trt1 -2.331441e-01  2.200422e-01
#> shape         1.783712e+00  3.031933e+00
#> Intercept     1.693988e+00  1.896088e+00
#> lprior       -4.862980e+00 -4.190513e+00
#> lp__         -6.679184e+02 -6.613081e+02
# }
```
